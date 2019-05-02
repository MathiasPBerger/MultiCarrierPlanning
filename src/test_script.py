from draft_data_class import Data
from draft_opti import model_builder
from draft_postprocess_class import PostProcess
from pyomo.opt import SolverFactory
from draft_print_class import Print
# from draft_plot_class import Plot
from os import makedirs
from os.path import isdir, join
from time import strftime

path = "../data/"
save_path = "../"

def init_folder():

    date = strftime("%Y%m%d")
    time = strftime("%H%M")
    
    if not isdir("../output"):
        makedirs("../output")
        path = "../output/run_"+date+'_'+str(time)
        makedirs(path)
    else:
        path = "../output/run_"+date+'_'+str(time)
        makedirs(path)
    
    return path
	
folder_path = init_folder()
		
data = Data(path)

model = model_builder(data, folder_path, write_lp=False)

opt = SolverFactory("gurobi")
# opt.options['MIPGap'] = 0.001
# #opt.options['Seed'] = 50000
# opt.options['Method'] = 2
# opt.options['OutputFlag'] = 1
# opt.options['LogFile'] = join(folder_path, 'gurobi.log')
#
results = opt.solve(model, tee=True, keepfiles=False)

printer = Print(model,data)
# # plotter = Plot(model, data, save_fig=False)
#
printer.print_output(folder_path)
# # plotter.plot(folder_path)
#
postprocessor = PostProcess(data, model)