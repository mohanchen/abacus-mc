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
source "$PROPS_SCRIPT_DIR/props_ml.sh"

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
# ML collectors (props_ml.sh): MLKEDF descriptors
#--------------------------------------------
run_ml_descriptor_props

#--------------------------------------------
# basic collectors that run after ml: imp_sol
#--------------------------------------------
run_basic_props_post_ml

#--------------------------------------------
# ML collectors (props_ml.sh): RPA
#--------------------------------------------
run_ml_rpa_props

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
# ML collectors (props_ml.sh): linear response excitations
#--------------------------------------------
run_ml_lr_props
#--------------------------------------------
# ML collectors (props_ml.sh): RDMFT energy terms
#--------------------------------------------
run_ml_rdmft_props

#--------------------------------------------
# basic collectors that run after rdmft: alllog
#--------------------------------------------
run_basic_props_post_rdmft

#--------------------------------------------
# trailing total-time entry
#--------------------------------------------
props_finalize "$1"
