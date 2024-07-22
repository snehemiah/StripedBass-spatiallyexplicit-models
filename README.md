# StripedBass-spatiallyexplicit-models

# Background
Updated folder with the final base version of the spatially-explicit, multi-stock models (.tpl) and data (.dat) for the Atlantic Striped Bass and varying models used for sensitivy runs developed by Samara Nehemiah, Ph.D. Student at the Chesapeake Biological Laboratory. All models were developed in ADMB (https://www.admb-project.org/). 

# Model Code
Code for the base model can be found in the **final_model** folder. This model is a two stock (Atlantic Coast and Chesapeake Bay), two area (Ocean and Chesapeake Bay) model that estimates abundance, movement, fishing mortality, and biomass at two 6-month time-steps from 1982-2017. This model assumes logisitc selectivitiy in the Ocean fishery and double logistic in the Bay fishery. There are 7 fisheries independey adult surveys used as data inputs, including 2 in the Bay and 5 in the Ocean. There are two age-1 surveys, including one in each region. There are 3 young-of-the-year surveys, including two in the Ocean region and one aggregate survey in the Bay. This model assumes aging bias in age composition data for the fisheries and certain surveys, as most agencies utilize scales to inform ages. 

Sensitivity runs and their assumptions can be found in the following folders:
1. **Final_time_vary_m** - this model assumes the same assumptions as the base model. However, this model allows natural mortality to vary overtime, based on the results of [Schonfeld (2023)](https://scholarworks.wm.edu/etd/1686662915/) which estimated natural mortality for Atlantic Striped Bass increasing every decade from 1990-2020.
2. **Final_without_agingerror** - this model assumes the same assumptions as the base model except it does not account for potential aging bias from scale-aged fish in age compoisition of fisheries dependent and independent data sources.
3. **Final_with_logistic_agingerror** - this model assumes the same assumptions as the base model, except is assumes logistic fishing selectivity in the Chesapeake Bay fishery.
4. **Final_with_logistic_noagingerror** - this model assumes the same assumptions as the base model, except is assumes logistic fishing selectivity in the Chesapeake Bay fishery and does not assume aging bias in age composition.
5. **Final_with_lowocc** - this model assumes the same assumptions as the base model, except it explore how changes to movement (occupancy probabilities) impact model estimates. 

# Partners
This work was funded by the NOAA Chesapeake Bay Office and the NMFS/Sea Grant Population and Ecosystems Dynamic Fellowship. Collaborators include researchers from the Chesapeake Biological Laboratory, Virginia Institute of Marine Science, and the Atlantic States Marine Fisheries Comission. Data for this research was provided by the Atlantic States Marine Fisheries Commission, Maine Department of Marine Resources, Rhode Island Department of Environmental Management, Massachusetts Division of Marine Fisheries, Connecticut Department of Energy and Environmental Protection, New York State Department of Environmental Conservation, New Jersey Division of Fish and Wildlife, Delaware Department of Natural Resources and Environmental Control, Maryland Department of Natural Resources, Potomac River Fisheries Commission, Virginia Institute of Marine Science, Virginia Department of Marine Resources, and North Carolina Department of Environmental Quality.

# Contact
If you have any questions, please contact Samara at snehemiah@umces.edu. 

<img src="https://www.umces.edu/sites/default/files/UMCES-CBL-logo.jpg" jsaction="" class="sFlh5c pT0Scc iPVvYb" style="max-width: 600px; height: 221px; margin: 0px; width: 557px;" alt="UMCES CBL logo.jpg | University of Maryland Center for Environmental Science" jsname="kn3ccd" aria-hidden="false">
