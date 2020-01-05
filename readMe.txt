There are 3 phenomena that we want to model. 
1.evaporation
The evaporation is modeled via mass transfer from water of the liquid phase to the air which only take places at the liquid-gas interface.
Therefore, we have a for-loop to loop over all cells of the liquid domain and determine if the cell VOF is smaller than some threshold (0.9). If so the cell is treated at the interface. If not, the mass transfer rate is set to 0. 

2.co2Loss due to interphase transfer
The CO2 interface transfer is modeled as linearized mass transfer. The the cell VOF is at [0.001 0.7], the cell is believed at the interface.

3.consumption by cells.
The consumption is modeled via the heterogeneous reaction:
1 HCO3(w)---->1H2o(w) or 0H2o(w)?
Since this phenomena onlu occurs at the liquid phase, we have a for loop to detrermine the cell VOF value. if the VOF is larger than some threshold (0.1), the cell is said to be at the liquid phase. Then the vertical location of the center of this cell is accessed via _CENTROID
and the light Depth is calculated. Moreover, the light intensity is calculated via the light intensity functiona and then, the reactionRate is interpolate from the lookup table. If the cell is not at the liquid phase, the reactionRate is set to be zero.

To achieve this, the UDF has the following functions:
1.DEFINE_EXECUTE_ON_LOADING(Initialize, libname) 
2.DEFINE_EXECUTE_AT_END(resetSpeciesFraction)
3.DEFINE_ON_DEMAND(resetSpeciesFractionByPH)
4.DEFINE_LINEARIZED_MASS_TRANSFER(co2LossV2, cell, thread, from_index,from_species_index, to_index, to_species_index, d_mdot_d_vof_from,d_mdot_d_vof_to)
5.DEFINE_HET_RXN_RATE(consumption, c, t, hr, mw, yi, rr, rr_t)
6. DEFINE_LINEARIZED_MASS_TRANSFER(evaporation, cell, thread, from_index,from_species_index, to_index, to_species_index, d_mdot_d_vof_from,d_mdot_d_vof_to)

I have changed the order of the species so that we can set the value of the proton and the mass fraction of water is calculated automatically.
The order is:
0:H   1
1:HCO3   63
2:CO3	62
3:CO2 	44
4:H2O 	18
