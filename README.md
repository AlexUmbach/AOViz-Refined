# What is AOViz?
AOViz is an R and RShiny based app for visualizing amplicon sequence data. It is designed to allow users to quickly and efficiently visualize their amplicon datasets and produce "near-publication-ready" figures. Importantly, AOViz provides users with a high level of interactivity and customization during figure generation. By seemlessly combining amplicon datasets with provided metadata, users can more readily explore and visualize microbial community compositional patterns. AOViz is, and will continue to be, a work in progress. The main branch is stable and should do everything it says it will, but if you encounter any errors please let me know and I'll do my best to address them in a timely manner. 

# Usage summary
Because AOViz is built using RShiny, users can simply download this repository to any location on their local machine and execute it by running the server.R script in RStudio (after installing the appropriate packages, of course).

Users require two files: an ASV table and an associated metadata table, both in a .txt or .tsv format. It is <b>critical</b> that the column names match the examples below <i>exactly</i>. Specifically, "Feature ID", "Consensus.Lineage", and "ReprSequence" in the ASV table, and "SampleName" in the metadata table. Note that the metadata "SampleName" column contains sample names that match the sample names in the ASV table. If you are having errors using AOViz, this is the first place to check. 

<b>ASV table</b>
![ASV_example](https://github.com/AlexUmbach/AOViz/assets/56092913/4af226dc-c0b0-4f96-a7bd-8ce8fac773f2)

<b>Metadata table</b>

![meta_example](https://github.com/AlexUmbach/AOViz/assets/56092913/92f7a906-9470-4a62-9d37-ee022e306062)

# Example figures
AOViz produces four main plots: a sequencing depth plot, a taxonomic bar plot, a relative abundance bubble plot, and a PCoA triplot. These four plots are contained within seprate RShiny tabs and run independent of each other. Descriptions of these plots are produces are contained in their respective tabs. For each tab, there is an option to download the figure and associated data tables, as well as statistics when applicable (i.e., triplot). The examples below are directly out of AOViz


# Read plot (bar and box)
![read_plot_bar_example](https://github.com/AlexUmbach/AOViz/assets/56092913/03e1096a-39ac-44a1-a8f2-b1eec417008b)
![read_plot_box_example](https://github.com/AlexUmbach/AOViz/assets/56092913/64021541-698e-452f-b9d0-767c1ca26089)

# Taxonomy bar plot (class)
![taxa_barplot_example](https://github.com/AlexUmbach/AOViz/assets/56092913/da04498a-6c86-47fc-83b0-d92caf03bbe9)

# Bubble plot
![bubbleplot_5ab_example](https://github.com/AlexUmbach/AOViz/assets/56092913/8cdd7aba-06f6-4865-860d-1e57c2b1ad68)

# Bray-Curtis triplot. 
For this plot, arrows represent numeric medata fit to the PCoA coordinates. These arrows are associated with p-values and R2 values (contained within a stats table available for download). Taxonomic relative abundance information is also mapped to PCoA data and shown in ordination space. Also available in UniFrac if provided a tree in newick format. 
![braycurtis_triplot_example](https://github.com/AlexUmbach/AOViz/assets/56092913/d0627b17-48ae-40dc-b6b6-41b96e659dfe)
