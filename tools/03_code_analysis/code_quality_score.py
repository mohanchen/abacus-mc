#!/usr/bin/env python3
"""ABACUS code quality scoring tool.

Scans source files/directories and assigns each file a quality score
starting from 100, deducting points for each rule violation found.

Rules (per-file score starts at 100):
  Core rules:
    - filename stem longer than 20 chars: -1
    - filename contains uppercase letters: -1
    - filename extension is .hpp: -50
    - each public member variable in a class/struct: -1
    - member function longer than 50 lines: -1 per additional 50-line block
  Zero-cost rules (no C++ parsing needed):
    - tab indentation: -1 per line (cap 5)
    - `using namespace std;`: -1 per occurrence (cap 5)
    - line longer than 120 chars: -1 per line (cap 5)
    - Chinese characters in comments/code: -1 per line (cap 5)
    - UPPERCASE constant naming (>3 chars all caps): -1 per occurrence (cap 5)
  Interface & dependency rules:
    - function declaration with default parameter: -2 per occurrence (cap 5)
    - GlobalV::/GlobalC::/PARAM.* cross-layer dependency: -3 per occurrence (cap 10)
    - #include of .hpp implementation header: -2 per occurrence (cap 5)
    - `friend` keyword exposing internals: -1 per occurrence (cap 5)
    - unpaired `new` without matching `delete`: -1 per occurrence (cap 5)
    - local variable shadowing a member variable: -1 per occurrence (cap 5)

Usage:
    python3 code_quality_score.py source/source_base
    # writes report to ./code_quality_score.txt (text default output file)
    python3 code_quality_score.py source/ -o report.txt
    python3 code_quality_score.py --format json path/to/file.cpp
    python3 code_quality_score.py --format json source/ -o report.json
    python3 code_quality_score.py --min-score 80 source/
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple


SOURCE_EXTENSIONS = {".c", ".cc", ".cpp", ".cxx", ".cu", ".h", ".hh", ".hpp", ".hxx", ".cuh"}

SKIP_DIRS = {
    ".git", "build", "__pycache__", "node_modules", ".cache",
    "third_party", "thirdparty", ".vscode", ".idea", ".trae-cn",
    "Dependencies",
}

CAPS = {
    "tab_indentation": 5,
    "using_namespace_std": 5,
    "line_too_long": 5,
    "chinese_comment": 5,
    "uppercase_constant": 5,
    "default_parameter": 5,
    "global_dependency": 10,
    "hpp_include": 5,
    "friend_keyword": 5,
    "unpaired_new_delete": 5,
    "member_local_name_conflict": 5,
}

WEIGHTS = {
    "filename_too_long": 1,
    "filename_uppercase": 1,
    "hpp_implementation": 50,
    "public_member_variable": 1,
    "member_function_too_long": 1,
    "tab_indentation": 1,
    "using_namespace_std": 1,
    "line_too_long": 1,
    "chinese_comment": 1,
    "uppercase_constant": 1,
    "default_parameter": 2,
    "global_dependency": 3,
    "hpp_include": 2,
    "friend_keyword": 1,
    "unpaired_new_delete": 1,
    "member_local_name_conflict": 1,
}

FUNCTION_LENGTH_THRESHOLD = 50
FUNCTION_LENGTH_STEP = 50
LINE_LENGTH_LIMIT = 120
FILENAME_LENGTH_LIMIT = 20
PASS_THRESHOLD = 60

CHINESE_RE = re.compile("[\u4e00-\u9fff]")
USING_NS_STD_RE = re.compile(r"\busing\s+namespace\s+std\b")
UPPERCASE_CONST_RE = re.compile(r"\b[A-Z][A-Z0-9_]{2,}\b")
ACCESS_RE = re.compile(r"^\s*(public|private|protected)\s*:")
CLASS_OPEN_RE = re.compile(r"\b(class|struct)\s+(\w+)\b")

GLOBAL_DEPENDENCY_RE = re.compile(r"\b(?:GlobalV::|GlobalC::|PARAM(?:\.|->|::))")
HPP_INCLUDE_RE = re.compile(r'^\s*#\s*include\s+[<"][^>"]+\.hpp[>"]')
FRIEND_RE = re.compile(r"\bfriend\b")
NEW_EXPR_RE = re.compile(r"\bnew\s+\w")
DELETE_EXPR_RE = re.compile(r"\bdelete\s*\[\s*\]?\s+\w")

DEFAULT_PARAM_RE = re.compile(
    r"\b\w+\s*\([^();]*\b\w+\s*=(?![=>])[^();]*\)\s*"
    r"(?:const\s*|noexcept\s*|override\s*|final\s*)*[;{]"
)

# Strict declaration regex: requires `type ... name ;` or `type ... name = ...;`.
# Excludes pure assignments like `counter = counter + 1;` which lack a leading type.
DECL_RE = re.compile(r"^\s*(?:\w+[\s\*]+)+(\w+)\s*(?:=[^;]*)?;\s*$")

BLOCK_KEYWORD_PREFIXES = (
    "public:", "private:", "protected:",
    "typedef", "using", "static_assert",
    "class", "struct", "friend",
    "//", "/*", "*",
    "return", "if ", "for ", "while ",
    "switch ", "case ", "template",
    "break", "continue", "goto", "default:",
    "else", "do", "try", "catch", "throw",
    "namespace", "enum", "union",
)


@dataclass
class Finding:
    rule: str
    line: Optional[int]
    reason: str
    deduction: int


@dataclass
class FileReport:
    path: str
    score: int
    findings: List[Finding] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "score": self.score,
            "findings": [
                {
                    "rule": f.rule,
                    "line": f.line,
                    "reason": f.reason,
                    "deduction": f.deduction,
                }
                for f in self.findings
            ],
        }


def discover_files(paths: Sequence[str]) -> List[Path]:
    """Walk input paths and return source files, skipping generated dirs."""
    result: List[Path] = []
    for p in paths:
        path = Path(p)
        if path.is_file():
            if path.suffix in SOURCE_EXTENSIONS:
                result.append(path)
        elif path.is_dir():
            for sub in path.rglob("*"):
                if not sub.is_file():
                    continue
                if any(part in SKIP_DIRS for part in sub.parts):
                    continue
                if sub.suffix in SOURCE_EXTENSIONS:
                    result.append(sub)
    return result


def strip_comments(content: str) -> str:
    """Return content with comments replaced by spaces, preserving line numbers.

    Handles // line comments and /* */ block comments. String literals are
    respected so that '/' inside strings is not mistaken for a comment.
    """
    out = []
    i = 0
    n = len(content)
    in_string = False
    in_char = False
    while i < n:
        c = content[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(content[i + 1])
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
            continue
        if in_char:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(content[i + 1])
                i += 2
                continue
            if c == "'":
                in_char = False
            i += 1
            continue
        if c == '"':
            in_string = True
            out.append(c)
            i += 1
            continue
        if c == "'":
            in_char = True
            out.append(c)
            i += 1
            continue
        if c == "/" and i + 1 < n:
            if content[i + 1] == "/":
                # line comment until newline
                j = content.find("\n", i)
                if j < 0:
                    j = n
                out.append(" " * (j - i))
                i = j
                continue
            if content[i + 1] == "*":
                # block comment until */
                j = content.find("*/", i + 2)
                if j < 0:
                    j = n
                else:
                    j += 2
                block = content[i:j]
                # preserve newlines
                out.append(re.sub(r"[^\n]", " ", block))
                i = j
                continue
        out.append(c)
        i += 1
    return "".join(out)


def find_class_blocks(code: str) -> List[Tuple[int, int, str, str]]:
    """Find top-level class/struct blocks. Returns list of
    (start_line_1indexed, end_line_1indexed, kind, name).

    Walks brace matching starting from each `class X {` / `struct X {` opener.
    """
    blocks: List[Tuple[int, int, str, str]] = []
    i = 0
    n = len(code)
    while i < n:
        m = CLASS_OPEN_RE.search(code, i)
        if not m:
            break
        # find the opening brace after the class/struct header
        # allow inheritance clauses: class X : public Y { ... }
        brace_pos = code.find("{", m.end())
        if brace_pos < 0:
            break
        # reject if there's a `;` before the brace (forward declaration)
        if ";" in code[m.end():brace_pos]:
            i = m.end()
            continue
        depth = 1
        j = brace_pos + 1
        in_string = False
        in_char = False
        while j < n and depth > 0:
            c = code[j]
            if in_string:
                if c == "\\":
                    j += 2
                    continue
                if c == '"':
                    in_string = False
                j += 1
                continue
            if in_char:
                if c == "\\":
                    j += 2
                    continue
                if c == "'":
                    in_char = False
                j += 1
                continue
            if c == '"':
                in_string = True
                j += 1
                continue
            if c == "'":
                in_char = True
                j += 1
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth == 0:
            start_line = code[:m.start()].count("\n") + 1
            end_line = code[:j].count("\n") + 1
            blocks.append((start_line, end_line, m.group(1), m.group(2)))
            i = j + 1
        else:
            i = m.end()
    return blocks


def detect_access_at_line(class_lines: List[str], target_idx: int, default_access: str) -> str:
    """Determine the active access specifier for line `target_idx` inside class."""
    access = default_access
    for k in range(target_idx + 1):
        m = ACCESS_RE.match(class_lines[k])
        if m:
            access = m.group(1)
    return access


def is_public_member_var(line: str) -> bool:
    """Heuristic: does this stripped code line look like a public member variable
    declaration (not a function, not a typedef, not a nested class)?"""
    s = line.strip()
    if not s.endswith(";"):
        return False
    if "(" in s or ")" in s:
        return False
    if "{" in s or "}" in s:
        return False
    if "=" in s and "(" in s:
        return False
    bad_prefixes = (
        "public:", "private:", "protected:",
        "typedef", "using", "static_assert",
        "class", "struct", "friend",
        "//", "/*", "*",
        "return", "if ", "for ", "while ", "switch ", "case ",
        "template",
    )
    if s.startswith(bad_prefixes):
        return False
    # skip macro-like lines (all caps with parens already excluded above)
    # require at least one identifier character
    if not re.search(r"[A-Za-z_]", s):
        return False
    return True


def _match_var_decl(stripped_line: str) -> Optional[str]:
    """If line looks like a variable declaration `type name;` or
    `type name = ...;`, return the variable name; else None.

    Excludes pure assignments (e.g. `counter = counter + 1;`) which lack a
    leading type token. This is a best-effort heuristic, not a full parser.
    """
    s = stripped_line.strip()
    if not s or s.startswith(BLOCK_KEYWORD_PREFIXES):
        return None
    m = DECL_RE.match(s)
    return m.group(1) if m else None


def analyze_class_blocks(content: str) -> Tuple[List[Finding], List[Finding], List[Finding]]:
    """Analyze class/struct blocks for public member variables, long member
    functions, and member/local name conflicts.

    Returns (public_member_findings, long_function_findings, name_conflict_findings).
    """
    pub_findings: List[Finding] = []
    long_func_findings: List[Finding] = []
    conflict_findings: List[Finding] = []

    stripped = strip_comments(content)
    lines = stripped.split("\n")

    for start_line, end_line, kind, name in find_class_blocks(stripped):
        default_access = "private" if kind == "class" else "public"
        body_start = start_line - 1
        body_end = end_line - 1
        class_lines = lines[body_start:body_end + 1]

        # pre-compute access per line (for public member variable rule)
        access_per_line: List[str] = []
        cur_access = default_access
        for k in range(len(class_lines)):
            m = ACCESS_RE.match(class_lines[k])
            if m:
                cur_access = m.group(1)
            access_per_line.append(cur_access)

        # pass 1: collect all member variable names (any access section).
        # depth starts at 0; the class header line brings it to 1. Only lines
        # that begin AND end at depth 1 are class-body member declarations
        # (lines that open a function body bring depth from 1 to 2).
        member_var_names: set = set()
        depth = 0
        for line in class_lines:
            prev_depth = depth
            depth += line.count("{") - line.count("}")
            if prev_depth == 1 and depth == 1:
                var_name = _match_var_decl(line)
                if var_name:
                    member_var_names.add(var_name)

        # pass 2: walk class body for findings
        depth = 0
        func_start_abs = -1
        for idx in range(len(class_lines)):
            line = class_lines[idx]
            abs_line = body_start + idx + 1  # 1-indexed absolute line
            prev_depth = depth

            # public member variable: only at depth 1 (class body, not in a function)
            if prev_depth == 1 and access_per_line[idx] == "public":
                if is_public_member_var(line):
                    pub_findings.append(Finding(
                        rule="public_member_variable",
                        line=abs_line,
                        reason=f"public member in {kind} {name}: {line.strip()}",
                        deduction=WEIGHTS["public_member_variable"],
                    ))

            # member/local name conflict: only at depth >= 2 (function body)
            if prev_depth >= 2:
                var_name = _match_var_decl(line)
                if var_name and var_name in member_var_names:
                    conflict_findings.append(Finding(
                        rule="member_local_name_conflict",
                        line=abs_line,
                        reason=(
                            f"local variable '{var_name}' in {kind} {name} "
                            f"shadows a member variable"
                        ),
                        deduction=WEIGHTS["member_local_name_conflict"],
                    ))

            # apply brace counting for this line
            depth = prev_depth + line.count("{") - line.count("}")

            # function body start: transition from depth 1 to >=2
            if prev_depth == 1 and depth >= 2 and func_start_abs < 0:
                func_start_abs = abs_line
            # function body end: transition from >=2 back to 1
            if prev_depth >= 2 and depth == 1 and func_start_abs > 0:
                func_end_abs = abs_line
                func_lines = func_end_abs - func_start_abs + 1
                if func_lines > FUNCTION_LENGTH_THRESHOLD:
                    excess = func_lines - FUNCTION_LENGTH_THRESHOLD
                    blocks = (excess + FUNCTION_LENGTH_STEP - 1) // FUNCTION_LENGTH_STEP
                    long_func_findings.append(Finding(
                        rule="member_function_too_long",
                        line=func_start_abs,
                        reason=(
                            f"member function in {kind} {name} spans {func_lines} lines "
                            f"(exceeds {FUNCTION_LENGTH_THRESHOLD})"
                        ),
                        deduction=blocks * WEIGHTS["member_function_too_long"],
                    ))
                func_start_abs = -1

    return pub_findings, long_func_findings, conflict_findings


def analyze_file(path: Path) -> FileReport:
    """Analyze one file and return its FileReport."""
    findings: List[Finding] = []
    name = path.name
    stem = path.stem
    suffix = path.suffix

    # filename rules
    if len(stem) > FILENAME_LENGTH_LIMIT:
        findings.append(Finding(
            rule="filename_too_long",
            line=None,
            reason=(
                f"filename '{name}' stem has {len(stem)} chars "
                f"(limit {FILENAME_LENGTH_LIMIT})"
            ),
            deduction=WEIGHTS["filename_too_long"],
        ))
    if any(c.isupper() for c in stem):
        findings.append(Finding(
            rule="filename_uppercase",
            line=None,
            reason=f"filename '{name}' contains uppercase letters",
            deduction=WEIGHTS["filename_uppercase"],
        ))
    if suffix == ".hpp":
        findings.append(Finding(
            rule="hpp_implementation",
            line=None,
            reason=".hpp implementation header prohibited (use .cpp + .h split)",
            deduction=WEIGHTS["hpp_implementation"],
        ))

    # read content
    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        findings.append(Finding(
            rule="read_error",
            line=None,
            reason=f"failed to read file: {e}",
            deduction=0,
        ))
        return FileReport(path=str(path), score=100, findings=findings)

    lines = content.split("\n")

    # line-based rules with caps
    tab_count = sum(1 for l in lines if "\t" in l)
    long_line_count = sum(1 for l in lines if len(l) > LINE_LENGTH_LIMIT)
    using_ns_count = sum(
        1 for l in lines if USING_NS_STD_RE.search(strip_comments(l))
    )
    chinese_count = sum(1 for l in lines if CHINESE_RE.search(l))

    stripped_content = strip_comments(content)
    upper_const_count = 0
    for l in stripped_content.split("\n"):
        upper_const_count += len(UPPERCASE_CONST_RE.findall(l))

    # interface & dependency rules
    global_dep_count = len(GLOBAL_DEPENDENCY_RE.findall(stripped_content))
    hpp_include_count = sum(1 for l in lines if HPP_INCLUDE_RE.match(l))
    friend_count = len(FRIEND_RE.findall(stripped_content))

    # default parameter: scan each stripped line for function-decl signature
    default_param_count = 0
    for l in stripped_content.split("\n"):
        default_param_count += len(DEFAULT_PARAM_RE.findall(l))

    # unpaired new/delete (file-level): rough heuristic, cap applied below
    new_count = len(NEW_EXPR_RE.findall(stripped_content))
    delete_count = len(DELETE_EXPR_RE.findall(stripped_content))
    unpaired_new = max(0, new_count - delete_count)

    def append_capped(rule: str, count: int) -> None:
        if count == 0:
            return
        cap = CAPS.get(rule)
        actual = min(count, cap) if cap else count
        cap_text = f" (capped at {cap})" if cap and count > cap else ""
        findings.append(Finding(
            rule=rule,
            line=None,
            reason=f"{count} occurrence(s){cap_text}",
            deduction=actual * WEIGHTS[rule],
        ))

    append_capped("tab_indentation", tab_count)
    append_capped("line_too_long", long_line_count)
    append_capped("using_namespace_std", using_ns_count)
    append_capped("chinese_comment", chinese_count)
    append_capped("uppercase_constant", upper_const_count)
    append_capped("default_parameter", default_param_count)
    append_capped("global_dependency", global_dep_count)
    append_capped("hpp_include", hpp_include_count)
    append_capped("friend_keyword", friend_count)
    append_capped("unpaired_new_delete", unpaired_new)

    # class-based rules
    pub_findings, long_func_findings, conflict_findings = analyze_class_blocks(content)
    findings.extend(pub_findings)
    findings.extend(long_func_findings)
    # name conflict findings already respect cap via per-file cap on the rule
    if len(conflict_findings) > CAPS.get("member_local_name_conflict", len(conflict_findings)):
        cap = CAPS["member_local_name_conflict"]
        conflict_findings = conflict_findings[:cap]
    findings.extend(conflict_findings)

    # compute score
    score = 100
    for f in findings:
        score -= f.deduction
    if score < 0:
        score = 0

    return FileReport(path=str(path), score=score, findings=findings)


def render_text(reports: List[FileReport], min_score: Optional[int]) -> str:
    """Render reports as plain text."""
    if not reports:
        return "No files to analyze.\n"

    visible = reports if min_score is None else [r for r in reports if r.score <= min_score]
    sorted_reports = sorted(visible, key=lambda r: r.score)

    out: List[str] = []
    width = 70
    out.append("=" * width)
    out.append("Code Quality Score Report")
    out.append("=" * width)
    out.append("")
    out.append(f"{'Score':>6}  {'File'}")
    out.append("-" * width)
    for r in sorted_reports:
        out.append(f"{r.score:>6}  {r.path}")
    out.append("")

    for r in sorted_reports:
        if not r.findings:
            continue
        out.append("-" * width)
        out.append(f"File: {r.path}  (score: {r.score})")
        out.append("-" * width)
        for f in r.findings:
            loc = f"line {f.line}" if f.line else "file"
            out.append(f"  [{f.rule}] {loc}: {f.reason} (-{f.deduction})")
        out.append("")

    total = len(reports)
    visible_n = len(visible)
    avg = sum(r.score for r in reports) / total if total else 0.0
    passing = sum(1 for r in reports if r.score >= PASS_THRESHOLD)
    out.append("=" * width)
    out.append(f"Files scanned:        {total}")
    out.append(f"Files shown:          {visible_n}")
    out.append(f"Average score:        {avg:.1f}")
    out.append(f"Passing (>= {PASS_THRESHOLD}):  {passing}/{total}")
    out.append("=" * width)
    return "\n".join(out) + "\n"


def render_json(reports: List[FileReport], min_score: Optional[int]) -> str:
    visible = reports if min_score is None else [r for r in reports if r.score <= min_score]
    return json.dumps({
        "summary": {
            "total_scanned": len(reports),
            "total_shown": len(visible),
            "average_score": (
                sum(r.score for r in reports) / len(reports) if reports else 0.0
            ),
            "passing": sum(1 for r in reports if r.score >= PASS_THRESHOLD),
            "pass_threshold": PASS_THRESHOLD,
        },
        "files": [r.to_dict() for r in sorted(visible, key=lambda r: r.score)],
    }, indent=2)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ABACUS code quality scoring tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("paths", nargs="+", help="file or directory paths to scan")
    parser.add_argument(
        "--format", choices=["text", "json"], default="text",
        help="output format (default: text)",
    )
    parser.add_argument(
        "--min-score", type=int, default=None,
        help="only show files with score <= this value in output",
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="write report to this file instead of stdout "
             "(default: code_quality_score.txt when format is text and "
             "no explicit output is given)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    files = discover_files(args.paths)
    if not files:
        print(f"No source files found in: {args.paths}", file=sys.stderr)
        return 1

    reports = [analyze_file(f) for f in files]

    if args.format == "json":
        rendered = render_json(reports, args.min_score)
    else:
        rendered = render_text(reports, args.min_score)

    out_path = args.output
    if out_path is None and args.format == "text":
        out_path = "code_quality_score.txt"

    if out_path:
        Path(out_path).write_text(rendered, encoding="utf-8")
        print(
            f"Report written to {out_path} "
            f"({len(reports)} files scanned, "
            f"{sum(1 for r in reports if r.score >= PASS_THRESHOLD)} passing)",
            file=sys.stderr,
        )
    else:
        print(rendered)
    return 0


if __name__ == "__main__":
    sys.exit(main())
