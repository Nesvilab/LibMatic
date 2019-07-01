##### RT alignment with a reference run
from common_funcs import (raise_if_not, str_to_path, unexpanduser_quote, list_as_shell_cmd, name_no_ext, strIII, os_fspath)
import itertools, os, subprocess, sys, fnmatch
import pandas as pd, numpy as np, pathlib
import pickle
import matplotlib
matplotlib.use('Agg')


if 0:
	sys.argv = ["%(prog)s", "./RT_align", "./iproph", "/data/dattam/PROJECTS/CoreFacility/PDLC/dda-lib-atcc-mm-R1/workdir/iproph","/data/teog/tpp5/bin/"]

has_DDA = sys.argv[3] != "none"

rtalign_data_directory = str_to_path(sys.argv[1])
dia_pepxml_directory = str_to_path(sys.argv[2])
dda_pepxml_directory = str_to_path(sys.argv[3]) if has_DDA else None
TPP_BIN = str_to_path(sys.argv[4])

rtalign_data_directory.mkdir(parents=True, exist_ok=False)

rt_dicts_file = rtalign_data_directory / "RT_dicts.pickle"
PEPTIDE_PROB = 0.9
abs_paths = [None if e is None else e.resolve() for e in [rtalign_data_directory, dia_pepxml_directory, dda_pepxml_directory, TPP_BIN / 'indexmzXML']]
print((TPP_BIN / 'indexmzXML').resolve(strict=True))
print("\n".join(str(e) for e in abs_paths))


CWD = os.getcwd()
os.chdir(os_fspath(rtalign_data_directory))


# dia_pepxml_directory = rtalign_data_directory / "iproph"
assert dia_pepxml_directory.exists()
pep_xml_rt_aligned_dir = rtalign_data_directory / "RT_aligned"




if has_DDA is True:
	assert dda_pepxml_directory.exists()
	dda_iproph_pep_xmls = sorted(fn.resolve() for fn in dda_pepxml_directory.glob("*.pep.xml"))
	assert len(dda_iproph_pep_xmls) > 0
	dda_pep_xml_rt_aligned_dir = rtalign_data_directory / "RT_aligned_DDA"

# rt_align_dir = data_directory / "RTalign"
iproph_pep_xmls = sorted(fn.resolve() for fn in dia_pepxml_directory.glob("*.pep.xml"))
assert len(iproph_pep_xmls) > 0
iproph_pep_xmls_fn_no_exts = [fn.name.split(".", 1)[0] for fn in iproph_pep_xmls]
endswith_Q123 = [fnmatch.fnmatchcase(f, "*_Q[123]") for f in iproph_pep_xmls_fn_no_exts]
[is_DIA_Umpire_output, ] = set(endswith_Q123)
iproph_pep_xmls_grouped = [(k, list(v)) for k, v in itertools.groupby(iproph_pep_xmls, lambda x: x.name.split(".")[0][:-len('_Q*')])] if \
	is_DIA_Umpire_output else \
	[(fn.name.split(".")[0], [fn]) for fn in iproph_pep_xmls]




## https://github.com/msproteomicstools/msproteomicstools/blob/master/msproteomicstoolslib/data_structures/modifications.py
## https://github.com/msproteomicstools/msproteomicstools/blob/master/analysis/spectral_libs/spectrast2tsv.py
import lxml.etree, csv, io, pathlib



def get_unimod_dict():
	## https://github.com/msproteomicstools/msproteomicstools/blob/master/msproteomicstoolslib/data_structures/modifications_default.tsv
	modifications_default_tsv = r"""modified-AA	TPP-nomenclature	Unimod-Accession	ProteinPilot-nomenclature	is_a_labeling	composition-dictionary
C	C[160]	4	[CAM]	FALSE	"{'H': 3 ,'C': 2 ,'N':1 , 'O': 1 }"
M	M[147]	35	[Oxi]	FALSE	{'O': 1 }
W	W[202]	35	[Oxi]	FALSE	{'O': 1 }
H	H[153]	35	[Oxi]	FALSE	{'O': 1 }
K	K[136]	259	[+08]	TRUE	"{'C' : -6 , '13C' : 6 , 'N' : -2 , '15N' : 2 }"
R	R[166]	267	[+10]	TRUE	"{'C' : -6 , '13C' : 6 , 'N' : -4 , '15N' : 4 }"
E	E[111]	27	[PGE]	FALSE	"{'H' : -2, 'O' : -1 }"
Q	Q[111]	28	[PGQ]	FALSE	"{'H' : -3, 'N' : -1 }"
C	C[143]	26	[PCm]	FALSE	"{'C' : 2, 'O' : 1}"
N-term	n[43]	5	[CRM]	FALSE	"{'C' : 1, 'H' : 1,'N' : 1,'O' : 1 }"
S	S[167]	21	[Pho]	FALSE	"{'H' : 1, 'O' : 3, 'P' : 1}"
T	T[181]	21	[Pho]	FALSE	"{'H' : 1, 'O' : 3, 'P' : 1}"
Y	Y[243]	21	[Pho]	FALSE	"{'H' : 1, 'O' : 3, 'P' : 1}"
N	N[115]	7	[Dea]	FALSE	"{'H' : -1, 'N': -1, ""O"": 1}"
Q	Q[129]	7	[Dea]	FALSE	"{'H' : -1, 'N': -1, ""O"": 1}"
C	C[149]	39	[XXX]	FALSE	"{'H' : 2, 'C' : 1, ""S"" : 1}"
D	D[131]	35	[Oxi]	FALSE	{'O': 1 }
K	K[144]	35	[Oxi]	FALSE	{'O': 1 }
Y	Y[179]	35	[Oxi]	FALSE	{'O': 1 }
R	R[172]	35	[Oxi]	FALSE	{'O': 1 }
N	N[130]	35	[Oxi]	FALSE	{'O': 1 }
P	P[113]	35	[Oxi]	FALSE	{'O': 1 }
C	C[119]	35	[Oxi]	FALSE	{'O': 1 }
N	N[317]	43	[XXX]	FALSE	"{'C': 8, 'H': 15, 'N' : 1, 'O' : 6 }"
N	N[349]	142	[XXX]	FALSE	"{'C': 14, 'H': 23, 'N' : 1, 'O' : 9 }"
""".rstrip()
# print(repr("\n".join(modifications_default_tsv.splitlines()[1:])))
	return {(l[0],int(l[1].split("[",1)[1][:-1])):int(l[2]) for l in csv.reader(io.StringIO("\n".join(modifications_default_tsv.splitlines()[1:])),dialect="excel-tab")}
to_unimod = get_unimod_dict()

def get_scannum_to_rt(msms_file: pathlib.Path):
	assert msms_file.suffix.casefold() == '.mzXML'.casefold()
	from xml.dom import pulldom
	doc = pulldom.parse(os_fspath(msms_file))
	scannums, rts = [], []
	for event, node in doc:
		if event is pulldom.START_ELEMENT and node.tagName == "scan":
			scannum = int(node.getAttribute("num"))
			rt_str = node.getAttribute("retentionTime")
			assert rt_str.startswith("PT") and rt_str.endswith("S")
			# scannum_to_rt.append((scannum, float(rt_str[2:-1]))
			scannums.append(scannum)
			rts.append(float(rt_str[2:-1]))
	doc.stream.close()

	scannum_to_rt = np.empty((max(scannums) + 1,), dtype=np.float32)
	scannum_to_rt.fill(np.nan)
	for scannum, rt in zip(scannums, rts):
		scannum_to_rt[scannum]=rt

	return scannum_to_rt

def get_pep(e, scannum_to_rt):
	retention_time_sec = e.get("retention_time_sec")
	start_scan, end_scan = int(e.get("start_scan")), int(e.get("end_scan"))
	assert start_scan==end_scan
	scan_number = start_scan
	PrecursorCharge = int(e.get("assumed_charge"))
	rt = None if retention_time_sec is None else float(retention_time_sec)

	[ee] = e.findall("{*}search_result/{*}search_hit")
	seq = ee.get("peptide")

	[interprophet_result] = ee.findall("{*}analysis_result/{*}interprophet_result") # TODO should work on all pep.xmls
	probablity = float(interprophet_result.get("probability"))

	modification_info_p = ee.findall("{*}modification_info")
	if modification_info_p != []:
		[modification_info] = modification_info_p
		nterm_mod = modification_info.get("mod_nterm_mass")
		mods = ([] if nterm_mod is None else [(0, float(nterm_mod))]) + \
			   [(int(mod.get("position")), float(mod.get("mass")))
				for mod in ee.findall("{*}modification_info/{*}mod_aminoacid_mass")]
	else:
		mods = []
	unimods = {pos: "(UniMod:{})".format(to_unimod.get(("N-term" if pos == 0 else seq[pos - 1], int(round(mass))), (seq[pos - 1], mass)))
			   for pos, mass in mods}
	FullUniModPeptideName=''.join(aa + unimods.get(idx, "")  for idx, aa in enumerate([""]+list(seq)+[""]))
	["FullUniModPeptideName", "PrecursorCharge", "Tr_recalibrated"]
	return (probablity, FullUniModPeptideName, PrecursorCharge, scannum_to_rt[scan_number])
	return (seq, rt, probablity, mods, FullUniModPeptideName)

import typing

def get_mzXMLs_from_pep_xml(pepxml_file: typing.TextIO):
	REC = re.compile(' base_name="(.+?)"')
	currpos = pepxml_file.tell()
	t = pepxml_file.read()
	pepxml_file.seek(currpos)
	paths = map(pathlib.Path, REC.findall(t))
	return [path.with_suffix(".mzXML") for path in paths if path.is_absolute()]


def table_from_pep_xml(infile: pathlib.Path):
	tree = lxml.etree.parse(os_fspath(infile))

	spectrum_paths = tree.findall("/{http://regis-web.systemsbiology.net/pepXML}msms_run_summary")
	try:
		(msms_file,) = set(pathlib.Path(spectrum_path.get("base_name")).with_suffix(".mzXML").resolve(strict=True) for spectrum_path in spectrum_paths)
	except FileNotFoundError as e:
		[spectrum_path] = tree.findall(
			# "/{http://regis-web.systemsbiology.net/pepXML}msms_pipeline_analysis"
			"/{http://regis-web.systemsbiology.net/pepXML}msms_run_summary"
			"/{http://regis-web.systemsbiology.net/pepXML}search_summary"
			"/{*}parameter[@name='spectrum, path']"
		)
		msms_file = pathlib.Path(spectrum_path.get("value")).resolve(strict=True)
		import re
		infile.write_text(
			re.compile('base_name="(.+?)" ').sub(f'base_name="{msms_file}" ', infile.read_text('utf-8'))
		)

	scannum_to_rt = get_scannum_to_rt(msms_file)
	gen = (get_pep(ee, scannum_to_rt) for ee in tree.findall("/{*}msms_run_summary/{*}spectrum_query"))
	p = set((FullUniModPeptideName, PrecursorCharge, rt) for probablity, FullUniModPeptideName, PrecursorCharge, rt in gen if probablity > PEPTIDE_PROB)
	colnames = ["FullUniModPeptideName", "PrecursorCharge", "Tr_recalibrated"]
	return pd.DataFrame({colname: e for colname, *e in zip(colnames, *p)})

###### get the RT alignment functions

def get_seq_charge_points(t):
	"""get dict of peptide sequence and unified charge for a tsv"""

	def unify_RTs(ch_rt_list):
		if len(ch_rt_list) == 1:
			return ch_rt_list[0][1]
		rts = [x[1] for x in ch_rt_list]
		# ptp = max(rts) - min(rts)
		# if ptp > 3:
		# 	return None
		# return np.mean(rts)
		med = np.median(rts)
		rts1 = [x for x in rts if abs(x - med) < 2]
		if len(rts1) == 0:
			return None
		return np.mean(rts1)

	p=set(map(tuple,t[["FullUniModPeptideName","PrecursorCharge","Tr_recalibrated"]].values))
	# fullname_charges=set(map(tuple,t[["FullUniModPeptideName","PrecursorCharge"]].values))
	# len(fullname_charges)
	# len(p)-len(fullname_charges)

	## combine peps with different charges
	d = dict()
	for fullname, pepcharge, rt in p:
		# d.setdefault(fullname, dict()).setdefault(pepcharge,set()).add(rt)
		d.setdefault(fullname, []).append((pepcharge, rt))
		# d.setdefault(fullname, set()).add(rt)

	g=((seq, unify_RTs(ch_rt)) for seq, ch_rt in d.items())
	return {seq:rt for seq,rt in g if rt is not None}



import concurrent.futures
with concurrent.futures.ThreadPoolExecutor() as exe,\
	concurrent.futures.ProcessPoolExecutor() as exe2:
	def f(g):
		groupname, files = g
		# l = [table_from_pep_xml(infile) for infile in files]
		# l = list(map(table_from_pep_xml, files))
		l = list(exe2.map(table_from_pep_xml, files))
		p = sorted({ee for e in l for ee in map(tuple, e.values)})
		del l
		colnames = ["FullUniModPeptideName", "PrecursorCharge", "Tr_recalibrated"]
		t=pd.DataFrame({colname: e for colname, *e in zip(colnames, *p)})
		del p
		return get_seq_charge_points(t)

	# ds = list(exe.map(f, iproph_pep_xmls_grouped))
	ds_iter = exe.map(f, iproph_pep_xmls_grouped)

	if has_DDA:
		dda_ds = list(exe.map(lambda file: f((None, [file])), dda_iproph_pep_xmls))
	ds = list(ds_iter)


with rt_dicts_file.open("wb") as f:
	pickle.dump(ds, f)

#####
with rt_dicts_file.open("rb") as f:
	ds = pickle.load(f)

###  write tables
import functools, os.path
all_identified_peptides = sorted(functools.reduce(set.union,(set(d.keys()) for d in ds)))
commonprefix_length = len("interact-")\
	if os.path.commonprefix([groupname for groupname, files in iproph_pep_xmls_grouped]).startswith("interact-")\
	else 0


table_RT_all_peptides = pd.DataFrame(
	{groupname[commonprefix_length:]:np.array([d.get(pep,np.nan) for pep in all_identified_peptides],dtype=np.float32)
		for (groupname, files), d in zip(iproph_pep_xmls_grouped, ds)},
	index=all_identified_peptides
)
table_RT_all_peptides.to_csv("RT_table_all_peptides.tsv", sep="\t")

np.argmax(list(map(len,ds)))
import functools
common_peptide_set = sorted(functools.reduce(set.intersection,(set(d.keys()) for d in ds)))
if len(common_peptide_set) < 11:
	import pandas as pd, numpy as np

	# t = pd.read_csv('file:///home/ci/tmp/RT_table_all_peptides.tsv', sep='\t', index_col=0)
	t = table_RT_all_peptides
	thresh = np.sort(t.notna().sum(1))[-11]
	t2 = t.dropna(thresh=t.shape[1] * .9)
	# t2=t.dropna(thresh=thresh)
	t2.to_csv("RT_table_common_peptides_relaxed.tsv", sep="\t")
	if 0:
		# NotImplementedError
		# https://stackoverflow.com/questions/33058590/pandas-dataframe-replacing-nan-with-row-average
		t3 = t2.fillna(t2.mean(1), axis=1)
	t3 = t2.T.fillna(t2.mean(1)).T
	t3.to_csv("RT_table_common_peptides_relaxed_fillna.tsv", sep="\t")
	# impute ds
	assert t3.shape[1] == len(ds)
	for d, (_, e) in zip(ds, t3.iteritems()):
		for k, v in e.iteritems():
			d[k] = v
	common_peptide_set = sorted(functools.reduce(set.intersection, (set(d.keys()) for d in ds)))

common_peptide_rts = np.transpose([[d[pep] for d in ds]
	for pep in common_peptide_set])
table_RT_common_peptides = pd.DataFrame(
	{groupname[commonprefix_length:]:np.array([d[pep] for pep in common_peptide_set],dtype=np.float32)
		for (groupname, files), d in zip(iproph_pep_xmls_grouped, ds)},
	index=common_peptide_set
)
table_RT_common_peptides.to_csv("RT_table_common_peptides.tsv", sep="\t")

### select reference run
## description:
'''
1. get the RTs of the common peptides and calculate the mean for each common peptide
2. the reference run has the minimum sum of squared difference against the mean RTs
'''

mean_rts = np.mean(common_peptide_rts, axis=0)
ref_run_idx = np.argmin(np.mean((common_peptide_rts - mean_rts)**2, axis=1)) # select central run

ref_run = ds[ref_run_idx]

def reg(x, y):
	import numpy as np
	if np.array_equal(x, y):
		return "ref run"

	import numpy as np, sklearn.isotonic

	import statsmodels.nonparametric.smoothers_lowess
	xy = statsmodels.nonparametric.smoothers_lowess.lowess(y, x, 0.1, 5)
	start, end = int(len(x)*0.1), int(len(x)*0.9)
	xy = xy[start:end]
	x = x[start:end]
	ir = sklearn.isotonic.IsotonicRegression()
	minx, maxx = min(x), max(x)
	yhat = ir.fit_transform(x, xy[:, 1])
	min_yhat = yhat[0]
	# when new_x<min(x)
	r0 = min_yhat / minx
	# when new_x > maxx - xr

	return ((minx, maxx), (0, yhat[-1]/maxx), r0, (x, yhat))


def predict(newx, regobj):
	if regobj == "ref run":
		return newx
	((minx, maxx), (intercept_, coef_), r0, (x, yhat)) = regobj
	if minx <= newx < maxx:
		idx = np.searchsorted(x, newx, 'right') - 1
		y0, y1 = yhat[idx:idx + 2]
		x0, x1 = x[idx:idx + 2]
		return y0 + (y1 - y0) * (newx-x0)/(x1-x0)
	if newx < minx:
		return r0 * newx
	assert maxx <= newx
	return intercept_ + coef_ * newx



def get_pts(d):
	g = ((rt, ref_run.get(seq, None))
		 for seq, rt in d.items())
	return np.array(sorted((x, y) for x, y in g if y is not None))

def get_RT_alignment_function(pts):
	return reg(pts[:,0], pts[:,1])

##### rewrite the iproph pepxmls


import re


def write_RT_aligned_pepxml(iproph_pep_xml, rt_aligned_pepxml, reg_obj):
	# reg_obj_max_idx = len(reg_obj) - 1
	# def repl(x):
	# 	idx = round(float(x.group()) * 10)
	# 	return str(reg_obj[min(idx,reg_obj_max_idx)])
	if reg_obj == "ref run":
		with rt_aligned_pepxml.open("wt") as newf, \
				iproph_pep_xml.open("rt") as origf:
			import shutil
			shutil.copyfileobj(origf, newf)

	def repl(x):
		return str(predict(float(x.group()), reg_obj))
		## spectrast will fail reading mzXML file if the replacement is not of the same length
		return format(predict(float(x.group()), reg_obj), ".6f")[:len(x.group())]
		return format(predict(float(x.group()), reg_obj), ".2f")

	import pathlib, re, filecmp
	def get_msms_from_pep_xml(p):
		t = p.read_text()
		paths = [pathlib.Path(p) for p in re.compile('<msms_run_summary base_name="(.+?)"').findall(t)]
		paths1 = [(p.parent / p.stem).with_suffix('.mzXML') for p in paths]
		# (msms_file,) = set(filter(pathlib.Path.exists, paths1))
		msms_files = set(filter(pathlib.Path.exists, paths1))
		msms_file = next(iter(msms_files))
		assert all([filecmp.cmp(f, msms_file) for f in msms_files]), msms_files
		return msms_file

	msms_file = get_msms_from_pep_xml(iproph_pep_xml)
	recomp = re.compile('(?<=retention_time_sec=")(.+?)(?=")')
	recomp_base_name = re.compile('base_name="(.+?)"')
	with rt_aligned_pepxml.open("wt") as newf, \
			iproph_pep_xml.open("rt") as origf:

		recomp2 = re.compile('(?<=retentionTime="PT)(.+?)(?=S")')
		new_msms_file = rt_aligned_pepxml.parent / msms_file.name
		with msms_file.open("rt") as msms_file_obj, \
				new_msms_file.open("wt") as new_msms_file_obj:
			total_repl = 0
			for line_1 in msms_file_obj:
				line_new, count = recomp2.subn(repl, line_1)
				assert count in (0, 1)
				total_repl += count
				new_msms_file_obj.write(line_new)
			assert total_repl > 0

		proc = subprocess.Popen([os_fspath(TPP_BIN / 'indexmzXML'), os_fspath(new_msms_file)], stdout=subprocess.DEVNULL)
		for line in origf:
			if '<msms_run_summary' in line:
				newline, count = recomp_base_name.subn(f'''base_name="{msms_file.with_suffix('')}"''', line, 1)
				assert count == 1
				newf.write(newline)
			else:
				newf.write(recomp.sub(repl, line))
	proc.wait()
	assert proc.returncode == 0, [proc.args, proc.returncode]

	new_msms_file.unlink()
	(new_msms_file.parent / (new_msms_file.name + ".new")).rename(new_msms_file)


import concurrent.futures, multiprocessing
pts_list = list(map(get_pts, ds))
if has_DDA:
	dda_pts_list = list(map(get_pts, dda_ds))
# rt_alignment_regobjs = list(map(get_RT_alignment_function, pts_list))
with concurrent.futures.ProcessPoolExecutor(min(multiprocessing.cpu_count(), len(pts_list))) as exe:
	rt_alignment_regobjs = list(exe.map(get_RT_alignment_function, pts_list))
	if has_DDA:
		dda_rt_alignment_regobjs = list(exe.map(get_RT_alignment_function, dda_pts_list))
# https://github.com/scikit-learn/scikit-learn/issues/7720

rt_alignment_regobjs_pts = list(zip(rt_alignment_regobjs, pts_list))

oldfile_newfile_reg_list = [(iproph_pep_xml, pep_xml_rt_aligned_dir / iproph_pep_xml.name, rt_alignment_regobj)
 for (groupname, pep_xml_group), rt_alignment_regobj in zip(iproph_pep_xmls_grouped, rt_alignment_regobjs)
 for iproph_pep_xml in pep_xml_group]
# %%time
pep_xml_rt_aligned_dir.mkdir(exist_ok=True)
with concurrent.futures.ProcessPoolExecutor(multiprocessing.cpu_count() * 2) as exe:
	futs=[exe.submit(write_RT_aligned_pepxml, iproph_pep_xml, RT_aligned_pep_xml, rt_alignment_regobj)
		  for iproph_pep_xml, RT_aligned_pep_xml, rt_alignment_regobj in oldfile_newfile_reg_list]
{f.result() for f in futs}

if has_DDA:
	dda_rt_alignment_regobjs_pts = list(zip(dda_rt_alignment_regobjs, dda_pts_list))
	dda_oldfile_newfile_reg_list = [(iproph_pep_xml, dda_pep_xml_rt_aligned_dir / iproph_pep_xml.name, rt_alignment_regobj)
								for iproph_pep_xml,rt_alignment_regobj
									in zip(dda_iproph_pep_xmls, dda_rt_alignment_regobjs)]
	# %%time
	dda_pep_xml_rt_aligned_dir.mkdir(exist_ok=True)
	with concurrent.futures.ProcessPoolExecutor(multiprocessing.cpu_count() * 2) as exe:
		futs = [exe.submit(write_RT_aligned_pepxml, iproph_pep_xml, RT_aligned_pep_xml, rt_alignment_regobj)
				for iproph_pep_xml, RT_aligned_pep_xml, rt_alignment_regobj in dda_oldfile_newfile_reg_list]
	{f.result() for f in futs}


def get_goodness_of_fit(reg_obj, pts):
	# linreg_func, pts, reg_obj = get_RT_alignment_function(d)
	yhat = np.array([predict(e, reg_obj) for e in pts[:, 0]])
	return np.sqrt(((yhat - pts[:, 1]) ** 2).mean())  # RMS error
	return abs(yhat - pts[:, 1]).mean()


if has_DDA:
	gofs = [get_goodness_of_fit(reg_obj, pts) / 60 for reg_obj, pts in dda_rt_alignment_regobjs_pts]
	t = '\n'.join("{}\t{}".format(xml_path.name, gof) for xml_path, gof in zip(dda_iproph_pep_xmls, gofs))
	pathlib.Path("goodness_of_fit_DDA.tsv").write_text("filename\tstandard dev (minutes)\n"+t)


def plot_RT_curves(rt_alignment_regobjs_pts, fn, filename_stems):
	##### plot the pairwise RTs and curve fit
	import matplotlib.pyplot as plt
	# fig=plt.figure()
	# import math
	# plotrows=math.ceil(np.sqrt(len(ds)))
	# plotcols=math.ceil(len(ds)/plotrows)
	# fig, ax = plt.subplots(plotrows,plotcols,sharex = True, sharey = True)
	import matplotlib.backends.backend_pdf
	with matplotlib.backends.backend_pdf.PdfPages(fn) as pdf:
		# for idx, d in enumerate(ds):
		# for regobjs_pts in rt_alignment_regobjs_pts:
		assert len(filename_stems)==len(rt_alignment_regobjs_pts), (filename_stems, rt_alignment_regobjs_pts)
		for pepxml_group, regobjs_pts in zip(filename_stems, rt_alignment_regobjs_pts):
			# linreg_func, pts, reg_obj = get_RT_alignment_function(d)
			reg_obj, pts = regobjs_pts
			print(reg_obj)

			# plt.hist(np.fmin(np.abs(linreg.intercept+linreg.slope*pts[:,0] - pts[:,1]),100))
			# plt.xlim(0,100)

			# plotobj=ax[idx // plotrows, idx % plotrows]
			plt.xlim((0,10000));plt.ylim((0,10000))
			plt.axes().set_aspect('equal', 'datalim')
			plotobj = plt
			# plotobj=fig.add_subplot(plotrows,plotcols,idx+1)
			plt.ylabel("reference run retention time (seconds)")
			plt.xlabel("run retention time (seconds)")
			plotobj.title(os.fspath(pepxml_group))
			plotobj.scatter(pts[:, 0], pts[:, 1], alpha=0.1,linewidth=0)
			# plotobj.plot(pts[:, 0], linreg_func(pts[:, 0]), 'red')
			# plotobj.plot(pts[:, 0], reg_obj[np.round(np.array(pts[:, 0]*10)).astype(int)], 'red')
			# plotobj.plot(pts[:, 0], [predict(e,reg_obj) for e in pts[:,0]], 'red')
			x = np.linspace(0,10000,num=1<<10)
			plotobj.plot(x, [predict(e,reg_obj) for e in x], 'red')
			# plotobj.savefig("f{}.pdf".format(idx))
			pdf.savefig()
			# plt.close(plotobj)
			plt.close()
		#	np.polyfit(pts[:,0],pts[:,1],deg=1)

plot_RT_curves(rt_alignment_regobjs_pts, 'pairwise_RT_alignments.pdf', [e[0] for e in iproph_pep_xmls_grouped])
if has_DDA:
	plot_RT_curves(dda_rt_alignment_regobjs_pts, 'pairwise_RT_alignments_dda.pdf', [e.name.split(".")[0] for e in dda_iproph_pep_xmls])

os.chdir(CWD)


'''
rewrite the precursorIntensity to circumvent the Spectrast bug.
ERROR PEPXML IMPORT: Precursor M/Z: 2.03492e+06 outside range [0, 1e+06]. Scan not imported.
PEPXML IMPORT: Problem loading spectrum. Skipped query "JHU_LM_DIA_Pancreatic_Pooled_DDA_F1_01.03030.03030".
'''

import pathlib, re, concurrent.futures, typing
mzXMLs = sorted(pathlib.Path().glob('RT_align/RT_aligned_DDA/*.mzXML'))
rec = re.compile('precursorIntensity="(.+?)"')
def repl(mo:typing.Match[str]):
	fill = '' if len(mo[1])==1 else '.'+'0'*(len(mo[1])-2)
	return f'precursorIntensity="0{fill}"'
def func(f:pathlib.Path):
	f2 = f.with_name(f.name+"0")
	with f.open() as fo, f2.open('w') as fo2:
		for line in fo:
			if rec.search(line):
				line = rec.sub(repl, line, 1)
				# line = rec.sub('precursorIntensity="0"', line, 1)
			fo2.write(line)
	if 1:
		f.rename(f.with_name(f.name+'.orig'))
		f2.rename(f)
	print(f)

with concurrent.futures.ProcessPoolExecutor() as pool:
	for f in mzXMLs:
		pool.submit(func, f)

