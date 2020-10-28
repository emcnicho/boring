import openmdao.api as om

from geometry import SizeComp
from fluid_properties import FluidPropertiesComp
from vapor_thermal_resistance import VapThermResComp

if __name__ == "__main__":
    p = om.Problem()
    model = p.model
    nn = 1


    p.model.add_subsystem(name = 'size',
                          subsys = SizeComp(num_nodes=nn),
                          promotes_inputs=['L_evap', 'L_cond', 'L_adiabatic', 't_w', 't_wk', 'D_od', 'D_v'],
                          promotes_outputs=['r_i', 'A_cond', 'A_evap', 'L_eff'])
    
    p.model.add_subsystem(name = 'fluids',
                          subsys = FluidPropertiesComp(num_nodes=nn),
                          promotes_inputs=['Q_hp', 'A_cond', 'h_c', 'T_coolant'],
                          # promotes_inputs=['A_cond'],
                          promotes_outputs=['T_cond', 'T_hp', 'T_hpfp', 'P_v', 'h_fg', 'rho_l','rho_v', 'mu_l', 'mu_v', 'k_l', 'k_v', 'sigma_l', 'cp_l', 'cp_v', 'v_fg', 'R_g'])
    
    p.model.add_subsystem(name = 'vapors',
                          subsys = VapThermResComp(num_nodes=nn),
                          promotes_inputs=['R_g', 'mu_v', 'T_hp', 'h_fg', 'P_v', 'rho_v', 'L_eff', 'D_v'],
                          promotes_outputs=['r_h', 'R_v'])
    

p.setup()

p.set_val('L_evap',0.01)
p.set_val('L_cond',0.02)
p.set_val('L_adiabatic',0.03)
p.set_val('t_w',0.0005)
p.set_val('t_wk',0.00069)
p.set_val('D_od', 0.006)
p.set_val('D_v',0.00362)
p.set_val('Q_hp',1)
p.set_val('h_c',1200)
p.set_val('T_coolant',293)


p.run_model()
p.check_partials(includes='fluids',method ='cs', compact_print=True,show_only_incorrect=False)
p.check_partials(includes='size',method ='cs', compact_print=True,show_only_incorrect=False)
print('Finished Successfully')

print('\n', '\n')
print('--------------Outputs---------------')
# print('The r_h Value is.......... ', p.get_val('r_h'))
# print('The R_v Value is.......... ', p.get_val('R_v'))
print('\n', '\n')