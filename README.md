# StripedBass-spatiallyexplicit-models

# Background
This project aimed to develop a spatially-explicit statistical catch-at-age model to estimate abundance and fishing mortality rates of Atlantic Striped Bass in the Chesapeake Bay and along the Atlantic coast. This repository has the final versions of the multi-stock, spatially-explicit model code (.tpl) and dat file (.dat) for  Atlantic Striped Bass that were developed by Samara Nehemiah for her Ph.D. studies at the Chesapeake Biological Laboratory. All models were developed in ADMB (https://www.admb-project.org/). 

# Model Description
Each model is a two spawning stock (Atlantic Coast and Chesapeake Bay), two area (Ocean and Chesapeake Bay) model that estimates abundance, movement, fishing mortality, and biomass at two 6-month time-steps from 1982-2017. This model assumes logisitc selectivitiy in the Ocean fishery and double logistic in the Bay fishery. There are 7 fisheries independey adult surveys used as data inputs, including 2 in the Bay and 5 in the Ocean. There are two age-1 surveys, including one in each region. There are 3 young-of-the-year surveys, including two in the Ocean region and one aggregate survey in the Chesapeake Bay. The four models evalauted different in their assumptions of natural mortality and ageing error. 

# Code Available
Code is available for four different version of the spatially-explicit model and can be found in the folder names indicated below. 
1. **Final_no_ageingerror** -The first model assumes constant natural mortality over time and regions, and does not correct for ageing bias from scale-aged fish in age compoisition of fishery dependent and independent data sources.
2. **Final_ageingerror** - The second model, assumes constant natural mortality over time and regions, and does correct for ageing bias from scale-aged fish in age compoisition of fishery dependent and independent data sources. 
3. **Final_timem_no_ageingerror** - The third model assumes natural mortality varies between regions and over time, based on the results of [Schonfeld (2023)](https://scholarworks.wm.edu/etd/1686662915/) which estimated natural mortality for Atlantic Striped Bass increasing every decade from 1990-2020. The third model does not correct for ageing bias from scale-aged fish in age compoisition of fishery dependent and independent data sources.
4. **Final_timem_ageingerror** - The fourth model assumes natural mortality varies between regions and over time, and does correct for ageing bias from scale-aged fish in age compoisition of fishery dependent and independent data sources.

# Partner Institutions
This work was funded by the NOAA Chesapeake Bay Office and the NMFS/Sea Grant Population and Ecosystems Dynamic Fellowship. Collaborators include researchers from the Chesapeake Biological Laboratory, Virginia Institute of Marine Science, and the Atlantic States Marine Fisheries Comission. Data for this research was provided by the Atlantic States Marine Fisheries Commission, Maine Department of Marine Resources, Rhode Island Department of Environmental Management, Massachusetts Division of Marine Fisheries, Connecticut Department of Energy and Environmental Protection, New York State Department of Environmental Conservation, New Jersey Division of Fish and Wildlife, Delaware Department of Natural Resources and Environmental Control, Maryland Department of Natural Resources, Potomac River Fisheries Commission, Virginia Institute of Marine Science, Virginia Department of Marine Resources, and North Carolina Department of Environmental Quality.

# Manuscript Citation
Nehemiah, S., R. Latour, A. Schonfeld, K. Drew, G. Nelson, D. Secor, and M.J. Wilberg. Effects of ageing error and natural mortality assumptions on spatially-explicit estimates of abundance and fishing mortality of striped bass. In preparation.

# Contact
If you have any questions, please contact Samara (she/her) at snehemiah@umces.edu. 

<img src="https://www.umces.edu/sites/default/files/UMCES-CBL-logo.jpg" jsaction="" class="sFlh5c pT0Scc iPVvYb" style="max-width: 600px; height: 221px; margin: 0px; width: 557px;" alt="UMCES CBL logo.jpg | University of Maryland Center for Environmental Science" jsname="kn3ccd" aria-hidden="false">
