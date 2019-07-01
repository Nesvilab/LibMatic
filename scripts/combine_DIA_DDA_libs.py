##### combining spectral libraries from DDA and DIA
from common_funcs import (raise_if_not, str_to_path, unexpanduser_quote, list_as_shell_cmd, name_no_ext, strIII, os_fspath)

import pandas as pd, numpy as np
import itertools
import pathlib
import sys


if not True:
	sys.argv = ["%(prog)s", "/data/dattam/PROJECTS/DIA-DataAcquisition-Yue/workdir/con_lib.tsv", "/data/dattam/PROJECTS/DDA-DataAcquisition-Yue/con_lib.tsv"]
	sys.argv = ["%(prog)s",
				"/data/dattam/PROJECTS/CoreFacility/PDLC/dia-atcc-mm-R1/workdir/con_lib.tsv",
				# "/data/dattam/PROJECTS/CoreFacility/PDLC/dia-atcc-R1/workdir/libgen/con_lib.tsv",
				"/data/dattam/PROJECTS/CoreFacility/PDLC/dda-lib-atcc-mm-R1/workdir/con_lib.tsv"
				# "/data/dattam/PROJECTS/CoreFacility/PDLC/dda-lib-atcc-mm-singleInjection-R1/workdir/con_lib.tsv"
				# "/data/teog/LFQbench/TTOF6600_SWATH_1ug_64windows_subset/workdir_tpp5/libgen/con_lib.tsv",
				# "/data/teog/LFQbench/TTOF6600_SWATH_1ug_64windows_subset/workdir_tpp5/libgen_DDA/con_lib.tsv"
				# "/data/teog/LFQbench/TTOF6600_SWATH_1ug_32windows/workdir/con_lib.tsv",
				# "/data/teog/LFQbench/DDA/workdir/lib_TTOF6600_32win/con_lib.tsv",

				# "/data/teog/LFQbench/TTOF6600_SWATH_1ug_64windows/workdir/con_lib.tsv",
				# "/data/teog/LFQbench/DDA/workdir/lib_TTOF6600_64win/con_lib.tsv",

				# "/data/dattam/PROJECTS/CoreFacility/PDLC/dia-atcc-mm-R1/workdir/con_lib.tsv",
				# "/data/dattam/PROJECTS/CoreFacility/PDLC/dda-lib-atcc-mm-R1/workdir/con_lib.tsv"

				#"/data/dattam/PROJECTS/CoreFacility/PDLC/dda-lib-atcc-mm-singleInjection-R1/workdir/con_lib.tsv"
				]

dia_con_lib_path = str_to_path(sys.argv[1])
dda_con_lib_path = str_to_path(sys.argv[2])

assert dia_con_lib_path.exists()
assert dda_con_lib_path.exists()

dia = pd.read_csv(dia_con_lib_path, sep='\t')
dia_transition_group_id = dia["transition_group_id"].values
# last int of the trans group id
dia_trans_grp_id_end = int(dia_transition_group_id[-1].split("_")[0])
dia_startidx = np.fromiter(
	(
		next(v)
		for k, v in itertools.groupby(dia.index, lambda x: dia_transition_group_id[x])),
	dtype=dia.index.dtype)
dia_pep_name_charge = dia[["FullUniModPeptideName", "PrecursorCharge"]].values
dia_pep_name_charge_set = {tuple(dia_pep_name_charge[lo]) for lo in dia_startidx}
assert set(map(tuple, dia_pep_name_charge)) == dia_pep_name_charge_set

dda = pd.read_csv(dda_con_lib_path, sep='\t')
dda_pep_name_charge = dda[["FullUniModPeptideName", "PrecursorCharge"]].values
dda_pep_name_charge_set = set(map(tuple, dda_pep_name_charge))
only_dda = dda_pep_name_charge_set - dia_pep_name_charge_set
dda_transition_group_id = dda["transition_group_id"].values
dda_trans_grp_id_int = np.array([int(e.split("_", 1)[0]) for e in dda_transition_group_id]) + dia_trans_grp_id_end
# dda["transition_group_id"] = np.core.defchararray.add(
# 	np.core.defchararray.add(dda_trans_grp_id_int.astype(str), "_"),
# 	[e.split("_", 1)[1] for e in dda_transition_group_id])
dda_startidxs = np.fromiter(
	(
		next(v)
		for k, v in itertools.groupby(dda.index, lambda x: dda_transition_group_id[x])),
	dtype=dda.index.dtype)
dda_endidxs = np.append(dda_startidxs[1:], dda.shape[0])
dda_pep_name_charge_dict = {tuple(dda_pep_name_charge[lo]): (lo, up) for lo, up in zip(dda_startidxs, dda_endidxs)}
add_dda = sorted(dda_pep_name_charge_dict[pep_name_charge] for pep_name_charge in only_dda)
## add the dda df
dda_to_add = dda.iloc[np.concatenate([np.r_[lo:up] for lo, up in add_dda])]

combined = pd.concat([dia, dda_to_add])
combined.to_csv("con_lib_DIA_DDA_combined.tsv", sep="\t", index=False)

assert set(map(tuple,combined[["FullUniModPeptideName", "PrecursorCharge"]].values))==(dda_pep_name_charge_set | dia_pep_name_charge_set)

# print("union", len(dda_pep_name_charge_set | dia_pep_name_charge_set))
# print("intersection", len(dda_pep_name_charge_set & dia_pep_name_charge_set))
# print("set difference DDA-DIA", len(dda_pep_name_charge_set - dia_pep_name_charge_set))
# print("set difference DIA-DDA", len(dia_pep_name_charge_set - dda_pep_name_charge_set))


import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
figure, axes = plt.subplots(2)
import matplotlib_venn
import matplotlib.backends.backend_pdf
with matplotlib.backends.backend_pdf.PdfPages("con_lib_venn.pdf") as pp:
	axes[0].set_title("Peptides in spectral libraries")
	matplotlib_venn.venn2([dda_pep_name_charge_set, dia_pep_name_charge_set], set_labels = ("DDA", "DIA"), ax=axes[0])

	axes[1].set_title("Proteins (Razor peptides) in spectral libraries")
	# matplotlib_venn.venn2([set(dda["razor_Protein"]), set(dia["razor_Protein"])], set_labels = ("DDA", "DIA"), ax=axes[1])
	matplotlib_venn.venn2([set(dda["Protein"]), set(dia["Protein"])], set_labels = ("DDA", "DIA"), ax=axes[1])
	# plt.show()
	pp.savefig()
