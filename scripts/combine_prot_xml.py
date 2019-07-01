from common_funcs import (raise_if_not, str_to_path, unexpanduser_quote, list_as_shell_cmd, name_no_ext, strIII, os_fspath)
import itertools, os, subprocess, sys, fnmatch
import pandas as pd, numpy as np, pathlib
import pickle

if not True:
	sys.argv = ["%(prog)s", "./workdir/libgen/combined_prots", "./workdir/iproph", "/data/dattam/PROJECTS/CoreFacility/PDLC/dda-lib-atcc-mm-R1/workdir/iproph", "/data/teog/tpp5/bin/"]
	sys.argv = ["%(prog)s", "./DIA/workdir/libgen/combined_prots", "./DIA/workdir/iproph", "./DDA/workdir/iproph", "/data/teog/tpp5/bin/"]

has_DDA = sys.argv[3] != "none"

combined_prot_data_directory = str_to_path(sys.argv[1])
dia_pepxml_directory = str_to_path(sys.argv[2])
dda_pepxml_directory = str_to_path(sys.argv[3]) if has_DDA else None
TPP_BIN = str_to_path(sys.argv[4])

dia_pep_xmls = list(map(os_fspath, dia_pepxml_directory.glob("*.iproph.pep.xml")))
dda_pep_xmls = list(map(os_fspath, dda_pepxml_directory.glob("*.iproph.pep.xml")))

raise_if_not(dia_pepxml_directory.exists(), "nonexistant DIA pep xml directory")
raise_if_not(dda_pepxml_directory.exists(), "nonexistant DDA pep xml directory")
raise_if_not(len(dia_pep_xmls) > 0, "no DIA pep xml found")
raise_if_not(len(dda_pep_xmls) > 0, "no DDA pep xml found")

combined_prot_data_directory.mkdir(exist_ok=True)

subprocess.run([os_fspath(TPP_BIN / "ProteinProphet")] +
			   dia_pep_xmls +
			   dda_pep_xmls +
			   [os_fspath(combined_prot_data_directory / "interact.prot.xml")] +
			   ["IPROPHET", "MINPROB0.9"],
			   cwd=combined_prot_data_directory)
