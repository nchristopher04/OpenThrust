#include <blowdownModel.h>
#include <thermo_functions.h>

double tankProps(double oxyMass, double Pc, double &T_Kelvin, double &TankPressure) {


	TankPressure = nox_vp(Temp);
	return 0;
}

/* Equilibrium (instantaneous boiling) tank blowdown model */
/* Empty tank of liquid nitrous */
void Nitrous_tank_liquid(void)
{
	double specVolume;
	double Chamber_press_bar_abs;
	double delta_outflow_mass, deltaQ, deltaTemp;
	double Enth_of_vap;
	double Spec_heat_cap;
	double tc;
	static double vaporized_mass = 0.0;
	/* blowdown simulation using nitrous oxide property calcs subroutines */
	Omdot_tank_outflow = mdot_tank_outflow;
	Enth_of_vap = nox_enthV(hybrid.tank_fluid_temperature_K); /* Get enthalpy (latent heat) of vaporisation */
	Spec_heat_cap = nox_CpL(hybrid.tank_fluid_temperature_K); /* Get specific heat capacity of the liquid nitrous */
															  /* Calculate the heat removed from the liquid nitrous during its vaporisation */
	deltaQ = vaporised_mass_old * Enth_of_vap;
	/* temperature drop of the remaining liquid nitrous due to losing this heat */
	deltaTemp = -(deltaQ / (hybrid.tank_liquid_mass * Spec_heat_cap));
	hybrid.tank_fluid_temperature_K += deltaTemp; /* update fluid temperature */
												  /* reality checks */
	if (hybrid.tank_fluid_temperature_K < (183))
	{
		hybrid.tank_fluid_temperature_K = (183); /* lower limit */
		hybrid.hybrid_fault = 1;
	}
	else if (hybrid.tank_fluid_temperature_K >309.15) //upper limit 36C
	{
		hybrid.tank_fluid_temperature_K = 309.15; 
		hybrid.hybrid_fault = 2;
	}
	/* get current nitrous properties */
	hybrid.tank_liquid_density = nox_Lrho(hybrid.tank_fluid_temperature_K);
	hybrid.tank_vapour_density = nox_Vrho(hybrid.tank_fluid_temperature_K);
	Chamber_press_bar_abs = hybrid.chamber_pressure_bar; /* Bar Abs */
														 /* calculate injector pressure drop and mass flowrate */
	mdot_tank_outflow = //
	
		/* integrate mass flowrate using Addams second order integration formula */
		/* Xn=X(n-1) + DT/2 * ((3 * Xdot(n-1) - Xdot(n-2)) */

		delta_outflow_mass = 0.5 * delta_time * (3.0 * mdot_tank_outflow - Omdot_tank_outflow);
	/* drain the tank based on flowrates only */
	hybrid.tank_propellant_contents_mass -= delta_outflow_mass; /* update mass within tank for next iteration */

	specVolume = (1.0 / rhoLiq) - (1.0 / rhoVap);
	liquidMass= (tankVolume - (TotalMass / rhoVap)) / specVolume;
	vaporMass = TotalMass - liquidMass;
	
	/* update for next iteration */
	specVolume = old_liquid_nox_mass - liquidMass;

	/* Add a 1st-order time lag (of 0.15 seconds) to aid numerical stability (this models the finite time required for boiling) */
	tc = delta_time / 0.15;
	vaporized_mass = tc * (bob - vaporized_mass) + vaporized_mass; // 1st-order lag
	vaporised_mass_old = vaporized_mass;
}