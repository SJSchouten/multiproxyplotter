# Multiproxyplotter
# consists of three modules which calculates statistics PCA, MDS, clustering, Corr-matrix, from a given set of data

# SNF_LakeWarm_Stats_Script provides the stats, note that clusters must be saved in the excel sheet manually.
# SNF_LakeWarm_StackedPlot_Script_withbar provides a vertical stack plot with options to choose from often used for depth plots.
# SNF_LakeWarm_StackedPlot_Script_horizontal provides a horizontal stack plot often used for age plots.
# proxies must be provided in the correct format and are accompanied by the correct meta-data

line 1: empty, 
line 2: proxy name, 
line 3: seccondary proxy name (or interpretation), 
line 4: units line 5: colors of the graph [rgb] triplet 0-1 scaled.
line 5: plot type;  bar graphs 'b', line 'l', Scatter 's', boxes (clusters) 'r', area 'a', and stairs 't'.
line 6: on-or-off toggle, i.e. 1 is on and will be plotted 0 is off and will not be considerred.

