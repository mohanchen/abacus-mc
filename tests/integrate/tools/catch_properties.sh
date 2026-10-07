#!/bin/bash

# mohan add 2025-05-03
# this compare script is used in different integrate tests
#
# Entry point: sources the shared library (tool paths, helper functions,
# one-time INPUT switch parsing) and then runs each property block below.
PROPS_SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$PROPS_SCRIPT_DIR/props_common.sh"
props_init "$1"
source "$PROPS_SCRIPT_DIR/props_basic.sh"
source "$PROPS_SCRIPT_DIR/props_mat.sh"
source "$PROPS_SCRIPT_DIR/props_cube.sh"

# Property collectors run in the same order as the original monolithic
# script so result files stay byte-identical; props_finalize() writes the
# trailing totaltimeref entry. Add a new category by dropping a
# props_<cat>.sh next to the others and calling run_<cat>_props here.
run_basic_props

#-------------------------------
# matrix collectors: out_dm1 (props_mat.sh)
#-------------------------------
run_mat_dm1_props

#-------------------------------
# cube collectors: out_pot/out_elf (props_cube.sh)
#-------------------------------
run_cube_pot_props

#-------------------------------
# matrix collectors (props_mat.sh): S/H(k)/XC/eband/H(R)/NPZ/r/T/dH
#-------------------------------
run_mat_props

#---------------------------------------
# cube collectors (props_cube.sh): chg/tau/LDOS/wfc real-space/PW/LCAO
#---------------------------------------
run_cube_props

#--------------------------------------------
# matrix collector: out_dm (props_mat.sh)
#--------------------------------------------
run_mat_dm_props

#--------------------------------------------
# cube collectors (props_cube.sh): mulliken/pchg/fingerprints/spinor
#--------------------------------------------
run_cube_tail_props

#--------------------------------------------
# ML gene data descriptors (.npy)
#--------------------------------------------
descriptor_dir="OUT.autotest/MLKEDF_Descriptors"
if [ -d "$descriptor_dir" ]; then
	python3 $COLLECT_NPY_MEANS "$descriptor_dir" >> "$1"
fi

#--------------------------------------------
# basic collectors that run after ml: imp_sol
#--------------------------------------------
run_basic_props_post_ml

#--------------------------------------------
# random phase approximation
#--------------------------------------------
if ! test -z "$run_rpa" && [ $run_rpa == 1 ]; then
	Etot_without_rpa=`grep Etot_without_rpa log.txt | awk 'BEGIN{FS=":"} {print $2}' `
	echo "Etot_without_rpa $Etot_without_rpa" >> $1
	rpa_outdir=$(get_input_key_value "rpa_outdir" "INPUT")
	if [ -z "$rpa_outdir" ]; then
		rpa_outdir="./OUT.librpa"
	fi
	rpa_outdir=${rpa_outdir%/}
	shopt -s nullglob
	rpa_ref_files=(refcoulomb_*.txt refCs_*.txt refshrink_sinvS_*.txt)
	if [ ${#rpa_ref_files[@]} -gt 0 ]; then
		IFS=$'\n' rpa_ref_files=($(printf '%s\n' "${rpa_ref_files[@]}" | LC_ALL=C sort))
		unset IFS
		for onref in "${rpa_ref_files[@]}"; do
			oncal_name=${onref#ref}
			oncal="$rpa_outdir/$oncal_name"
			compare_key="CompareRPA_$(sanitize_result_key "$oncal_name")_pass"
			record_compare_result "$1" "$compare_key" "$onref" "$oncal" 8 1
		done
	fi
	shopt -u nullglob
fi

#--------------------------------------------
# deepks
#--------------------------------------------
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
bash ${script_dir}/catch_deepks_properties.sh $1

#--------------------------------------------
# basic collectors that run after deepks: symmetry
#--------------------------------------------
run_basic_props_post_deepks

#--------------------------------------------
# check currents in rt-TDDFT 
#--------------------------------------------
if ! test -z "$out_current" && [ $out_current ]; then
	current1ref=current_tot.txt.ref
	current1cal=OUT.autotest/current_tot.txt
	python3 $COMPARE_SCRIPT $current1ref $current1cal 10
	echo "CompareCurrent_pass $?" >>$1
fi

#--------------------------------------------
# Check electric fields in rt-TDDFT
#--------------------------------------------
if ! test -z "$out_efield" && [ "$out_efield" == 1 ]; then
	efield_refs=(efield_*.txt.ref)
	if [ ! -e "${efield_refs[0]}" ]; then
		echo "CompareEfieldReference_pass 1" >>$1
	else
		for efield_ref in "${efield_refs[@]}"; do
			efield_name=${efield_ref%.ref}
			efield_key=$(sanitize_result_key "$efield_name")
			record_compare_result "$1" "Compare${efield_key}_pass" "$efield_ref" "OUT.autotest/$efield_name" 8
		done
	fi
fi

#--------------------------------------------
# Check vector potential in rt-TDDFT
#--------------------------------------------
if ! test -z "$out_vecpot" && [ "$out_vecpot" == 1 ]; then
	record_compare_result "$1" "CompareVectorPot_pass" "vector_pot.txt.ref" "OUT.autotest/vector_pot.txt" 8
fi

#--------------------------------------------
# Linear response function 
#--------------------------------------------
if [ $is_lr == 1 ]; then
	shopt -s nullglob
	lr_files=(OUT.autotest/trans_analysis_*_tda.dat)
	if [ ${#lr_files[@]} -gt 0 ]; then
		cat "${lr_files[@]}" | awk '/Excitation Energy/{p=1; next} p && /^[[:space:]]*[0-9]+[[:space:]]/{printf "excitationenergyref%d %.6f\n", ++n, $2} /Occupied orbital/{p=0}' >>$1
	fi
	shopt -u nullglob
fi
#--------------------------------------------
# Check RDMFT method 
#--------------------------------------------
if ! test -z "$rdmft" && [[ $rdmft == 1 ]]; then
	echo "" >>$1
	echo "The following energy units are in Rydberg:" >>$1

	E_TV_RDMFT=$(grep "E_TV_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_TV_RDMFT_ref $E_TV_RDMFT" >>$1

	E_hartree_RDMFT=$(grep "E_hartree_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_hartree_RDMFT_ref $E_hartree_RDMFT" >>$1

	Exc_cwp22_RDMFT=$(grep "Exc_cwp22_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "Exc_cwp22_RDMFT_ref $Exc_cwp22_RDMFT" >>$1

	E_Ewald=$(grep "E_Ewald" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_Ewald_ref $E_Ewald" >>$1

	E_entropy=$(grep "E_entropy(-TS)" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_entropy_ref $E_entropy" >>$1

	E_descf=$(grep "E_descf" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_descf_ref $E_descf" >>$1

	Etotal_RDMFT=$(grep "Etotal_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "Etotal_RDMFT_ref $Etotal_RDMFT" >>$1

	Exc_ksdft=$(grep "Exc_ksdft" "$running_path" | tail -1 | awk '{print $2}')
	echo "Exc_ksdft_ref $Exc_ksdft" >>$1

	E_exx_ksdft=$(grep "E_exx_ksdft" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_exx_ksdft_ref $E_exx_ksdft" >>$1

	echo "" >>$1
fi

#--------------------------------------------
# basic collectors that run after rdmft: alllog
#--------------------------------------------
run_basic_props_post_rdmft

#--------------------------------------------
# trailing total-time entry
#--------------------------------------------
props_finalize "$1"
