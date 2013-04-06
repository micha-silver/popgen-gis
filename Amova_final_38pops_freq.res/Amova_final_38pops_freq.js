
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/arlequin31/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (Amova_final_38pops_freq.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 28/12/10 at 20:08:35", "Amova_final_38pops_freq.htm#28_12_10at20_08_35"))
	aux1 = insFld(foldersTree, gFld("Run of 28/12/10 at 20:44:11", "Amova_final_38pops_freq.htm#28_12_10at20_44_11"))
	aux1 = insFld(foldersTree, gFld("Run of 29/12/10 at 13:15:31", "Amova_final_38pops_freq.htm#29_12_10at13_15_31"))
	aux1 = insFld(foldersTree, gFld("Run of 29/12/10 at 13:57:38", "Amova_final_38pops_freq.htm#29_12_10at13_57_38"))
	aux1 = insFld(foldersTree, gFld("Run of 29/12/10 at 14:00:05", "Amova_final_38pops_freq.htm#29_12_10at14_00_05"))
	insDoc(aux1, gLnk("R", "Settings", "Amova_final_38pops_freq.htm#29_12_10at14_00_05_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure", "Amova_final_38pops_freq.htm#29_12_10at14_00_05_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "Amova_final_38pops_freq.htm#29_12_10at14_00_05_pairw_diff"))
		insDoc(aux2, gLnk("R", "Exact tests", "Amova_final_38pops_freq.htm#29_12_10at14_00_05_exct_tests"))
