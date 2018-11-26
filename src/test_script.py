from draft_data_class import Data
from draft_opti import model_builder
from draft_postprocess_class import PostProcess
from pyomo.opt import SolverFactory
from draft_print_class import Print
from draft_plot_class import Plot
from os import makedirs
from os.path import isdir, join
from time import strftime

ID = input('Write user here: ')

if ID == 'dz':
    path = "/home/dcradu/model/fluxys/1node/data/"
    save_path = "/home/dcradu/model/fluxys/1node/"

elif ID == 'dl': 
    path = "/Users/david/py_workspace/fluxys/model/data/" 
    save_path = "/Users/david/py_workspace/fluxys/model/" 
	
elif ID == 'ml':
    path = "/Users/mathiasberger/Documents/Programming/Python/MultiCarrierPlanning/data/"
    save_path = "/Users/mathiasberger/Documents/Programming/Python/MultiCarrierPlanning/"

elif ID == 'mz':
    path = "/home/berger/..."
    save_path = "/home/berger/..."


def init_folder():

    date = strftime("%Y%m%d")
    time = strftime("%H%M")
    path = "../output/run_"+date+'_'+str(time)
    if not isdir(path):
        makedirs(path)
    return path
	
folder_path = init_folder()
		
data = Data(path)

model = model_builder(data, folder_path, write_lp=True)

opt = SolverFactory("gurobi", solver_io='python')
opt.options['MIPGap'] = 0.05
#opt.options['Seed'] = 50000
opt.options['Method'] = 2
opt.options['OutputFlag'] = 1
opt.options['LogFile'] = join(folder_path, 'gurobi.log')

results = opt.solve(model, tee=True, keepfiles=False)

printer = Print(model,data)
#plotter = Plot(model, data, save_fig=False)

printer.print_output(folder_path)
#plotter.plot(folder_path)

postprocessor = PostProcess(data, model)