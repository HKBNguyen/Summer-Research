import os
import pprint
import re

pp = pprint.PrettyPrinter(indent=4)

####################################################################################################################

MEM_HEADER="""
universe = vanilla
notification = never
should_transfer_files = yes
when_to_transfer_output = always
copy_to_spool = false
executable = ../limit.sh
    """

MEM_JOB="""
requirements = regexp("slot({use_slots})@pedigree-({use_clusters}).ics.uci.edu", Name)
initialdir = {init_dir}
output = {uai_file}_{option_strings_name}.txt
error  = err_{uai_file}_{option_strings_name}.txt
log    = log_{uai_file}_{option_strings_name}.txt
transfer_input_files = {program}, {uai_loc}, {vo_loc}, {evid_loc}
arguments = {time_limit_shell} {memory_limit_kb} ./ARP-v3.0.2 -nR 2000000000 -igEH 0 -a context -fUAI {uai_file} -fVO {vo_file} {option_strings_cmd} -fEV {evid_file} -t {time_limit} -exactZ {Z_value}
queue

    """

####################################################################################################################

program_name = 'ARP-v3.0.2'
haplo1_root = 'Experiments' # e.g. haplo1_root / park / promedas / or_chacin... /
clusters_available = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]

local_benchmark_root = os.path.join(os.getcwd(), 'Experiments')
domains = ['DBN', 'Pedigree', 'Promedas','Grids']

command_per_solver = {}
#command_per_solver['aobb_or'] = '--algorithm aobb --heuristic wmb-mm --verbose --positive --seed 12345678 --orsearch'
#AO, OR, abstraction levels
command_per_solver['EXP_AO_0_N'] = '-treeType AO -nContext 0 -iB 10 -nLevelsLimit 999999 -proper 0'
command_per_solver['EXP_AO_4_N'] = '-treeType AO -nContext 4 -iB 10 -nLevelsLimit 999999 -proper 0'
command_per_solver['EXP_AO_8_N'] = '-treeType AO -nContext 8 -iB 10 -nLevelsLimit 999999 -proper 0'
#command_per_solver['EXP_AO_0'] = '-treeType AO -nContext 0 -iB 10'
#command_per_solver['EXP_AO_1'] = '-treeType AO -nContext 1 -iB 10'
#command_per_solver['EXP_AO_2'] = '-treeType AO -nContext 2 -iB 10'
ibounds = [10]
####################################################################################################################
# dictionary that searches for exact solution by problem name; if no exact solution (problem is large), then default value is 0
def foo (problem):
    Z_Dict ={
        #pedigrees
        "pedigree1": -14.11, 
        "pedigree13": -31.18,
        "pedigree18": -78.14,
        "pedigree19": -59.02,
        "pedigree20": -29.63,
        "pedigree23": -38.56,
        "pedigree25": -115.78,
        "pedigree30": -83.73,
        "pedigree31": -69.70,
        "pedigree33": -54.28,
        "pedigree34": -64.23,
        "pedigree37": -116.58,
        "pedigree38": -54.25,
        "pedigree39": -102.20,
        "pedigree40": -87.88,
        "pedigree41": -76.04,
        "pedigree42": -30.76,
        "pedigree44": -63.47,
        "pedigree50": -22.88,
        "pedigree51": -77.27,
        "pedigree7": -64.82,
        "pedigree9": -78.52,
        "Family2Dominant.1.5loci": -17.22,
        "Family2Dominant.20.5loci": -11.46,
        "Family2Recessive.15.5loci": -15.21,
        #grids
        "grid10x10.f10": 303.09,
        "grid10x10.f10.wrap": 333.32,
        "grid10x10.f15.wrap": 497.76,
        "grid10x10.f5.wrap": 169.41,
        "grid20x20.f10": 1311.98,
        "grid20x20.f15": 1962.98,
        "grid20x20.f2": 291.73,
        "grid20x20.f5": 665.12,
        #promedas
        "or_chain_1.fg": -10.76,
        "or_chain_10.fg": -8.39,
        "or_chain_102.fg": -10.04,
        "or_chain_106.fg": -10.06,
        "or_chain_107.fg": -11.57,
        "or_chain_111.fg": -5.86,
        "or_chain_12.fg": -12.22,
        "or_chain_128.fg": -10.91,
        "or_chain_129.fg": -9.42,
        "or_chain_132.fg": -8.43,
        "or_chain_133.fg": -6.71,
        "or_chain_138.fg": -10.05,
        "or_chain_139.fg": -7.00,
        "or_chain_140.fg": -5.54,
        "or_chain_147.fg": -8.14,
        "or_chain_148.fg": -8.13,
        "or_chain_149.fg": -8.60,
        "or_chain_15.fg": -11.85,
        "or_chain_150.fg": -7.25,
        "or_chain_152.fg": -10.15,
        "or_chain_153.fg": -8.69,
        "or_chain_154.fg": -1.19,
        "or_chain_155.fg": -9.95,
        "or_chain_161.fg": -9.24,
        "or_chain_165.fg": -21.80,
        "or_chain_168.fg": -1.80,
        "or_chain_17.fg": -3.16,
        "or_chain_175.fg": -1.27,
        "or_chain_176.fg": -2.20,
        "or_chain_18.fg": -4.05,
        "or_chain_180.fg": -2.81,
        "or_chain_182.fg": -3.08,
        "or_chain_186.fg": -9.47,
        "or_chain_188.fg": -12.31,
        "or_chain_197.fg": -1.19,
        "or_chain_198.fg": -4.30,
        "or_chain_201.fg": -1.74,
        "or_chain_203.fg": -1.25,
        "or_chain_209.fg": -11.46,
        "or_chain_214.fg": -14.97,
        "or_chain_218.fg": -1.49,
        "or_chain_220.fg": -10.33,
        "or_chain_225.fg": -1.74,
        "or_chain_236.fg": -4.17,
        "or_chain_242.fg": -5.35,
        "or_chain_247.fg": -4.98,
        "or_chain_248.fg": -11.55,
        "or_chain_29.fg": -11.93,
        "or_chain_38.fg": -7.51,
        "or_chain_4.fg": -11.72,
        "or_chain_42.fg": -3.34,
        "or_chain_45.fg": -6.98,
        "or_chain_48.fg": -4.79,
        "or_chain_53.fg": -9.18,
        "or_chain_61.fg": -8.00,
        "or_chain_62.fg": -9.46,
        "or_chain_64.fg": -8.84,
        "or_chain_65.fg": -4.67,
        "or_chain_68.fg": -4.34,
        "or_chain_72.fg": -7.06,
        "or_chain_74.fg": -5.58,
        "or_chain_85.fg": -2.49,
        "or_chain_90.fg": -8.55,
        "or_chain_91.fg": -11.32,
        "or_chain_93.fg": -8.72,
        #DBN 
        "rbm_20": 58.53,
        "rbm_21": 63.15,
        "rbm_22": 66.55,
        "rbm_ferro_20": 151.16,
        "rbm_ferro_21": 152.62,
        "rbm_ferro_22": 166.11,
        "rus_20_40_0_1": 617.31,
        "rus_20_40_0_2": 791.57,
        "rus_20_40_0_3": 903.04,
        "rus_20_40_1_1": 1004.68,
        "rus_20_40_1_2": 910.22,
        "rus_20_40_1_3": 900.64,
        "rus_20_40_2_1": 797.82,
        "rus_20_40_2_2": 689.28,
        "rus_20_40_2_3": 749.30,
        "rus_20_40_3_1": 825.86,
        "rus_20_40_3_2": 839.01,
        "rus_20_40_3_3": 841.92,
        "rus_20_40_4_1": 931.45,
        "rus_20_40_4_2": 935.72,
        "rus_20_40_4_3": 853.26,
        "rus_20_40_5_1": 875.91,
        "rus_20_40_5_2": 839.66,
        "rus_20_40_5_3": 881.95,
        "rus_20_40_6_1": 723.87,
        "rus_20_40_6_2": 982.42,
        "rus_20_40_6_3": 966.71,
        "rus_20_40_7_1": 752.12,
        "rus_20_40_7_2": 756.19,
        "rus_20_40_7_3": 869.75,
        "rus_20_40_8_1": 726.88,
        "rus_20_40_8_2": 898.68,
        "rus_20_40_8_3": 703.83,
        "rus_20_40_9_1": 866.61,
        "rus_20_40_9_2": 834.53,
        "rus_20_40_9_3": 864.01,
        #break
        "rus2_20_40_0_1": 106.92,
        "rus2_20_40_0_2": 120.62,
        "rus2_20_40_0_3": 88.20,
        "rus2_20_40_1_1": 107.43,
        "rus2_20_40_1_2": 104.48,
        "rus2_20_40_1_3": 101.62,
        "rus2_20_40_2_1": 119.26,
        "rus2_20_40_2_2": 118.57,
        "rus2_20_40_2_3": 112.44,
        "rus2_20_40_3_1": 110.14,
        "rus2_20_40_3_2": 141.20,
        "rus2_20_40_3_3": 123.54,
        "rus2_20_40_4_1": 99.16,
        "rus2_20_40_4_2": 107.05,
        "rus2_20_40_4_3": 95.18,
        "rus2_20_40_5_1": 95.46,
        "rus2_20_40_5_2": 88.32,
        "rus2_20_40_5_3": 169.89,
        "rus2_20_40_6_1": 119.14,
        "rus2_20_40_6_2": 68.35,
        "rus2_20_40_6_3": 134.12,
        "rus2_20_40_7_1": 88.48,
        "rus2_20_40_7_2": 111.55,
        "rus2_20_40_7_3": 81.44,
        "rus2_20_40_8_1": 106.90,
        "rus2_20_40_8_2": 99.86,
        "rus2_20_40_8_3": 100.59,
        "rus2_20_40_9_1": 94.55,
        "rus2_20_40_9_2": 94.58,
        "rus2_20_40_9_3": 63.61
    }
    return Z_Dict.get(problem, 0)
    #set memory limit as 8 GB, will test to see if more is needed
reg = re.compile('^Results*')
def run(solver_name, time_limit=3600, memory_limit=8000):
    instance_names = {}
    for domain in domains:
        domain_path = os.path.join(local_benchmark_root, domain)
        instance_names[domain] = sorted([f for f in os.listdir(domain_path) if os.path.isdir(os.path.join(domain_path, f)) and not reg.match(f)])
      #  results_folder = [f for f in os.listdir(domain_path) if reg.match(f)]

    ####################################################################################################################

    num_jobs_per_cluster = int(24000 / memory_limit)
    if num_jobs_per_cluster < 10:
        slot_assigned = '[1-' + str(num_jobs_per_cluster) + ']'
    elif num_jobs_per_cluster < 20:
        num_jobs = num_jobs_per_cluster - 10
        slot_assigned = '[1-9]|1[0-' + str(num_jobs) + ']'
    elif num_jobs_per_cluster < 30:
        num_jobs = num_jobs_per_cluster - 20
        slot_assigned = '[1-9]|1[0-9]|2[0-' + str(num_jobs) + ']'

    cluster_assigned = '|'.join( [str(c) for c in clusters_available] )

    ####################################################################################################################

    for domain in domains:
        domain_path = os.path.join(local_benchmark_root, domain)
        condor_script_name = os.path.join(domain_path, solver_name +'_' + domain +'.condor')
        condor_file = open(condor_script_name, 'w')
        condor_file.write(MEM_HEADER)

        for instance in instance_names[domain]:
            # instance_path = os.path.join(domain_path, instance)
            uai_file = instance + '.uai'
            evid_file = instance + '.uai.evid'
            vo_file = instance + '.uai.ord.elim'
    #        Z_Dict = {
     #           'pedigree1': -14.11,
      #          'pedigree13': -13.15

       #     }            
        #    Exact_Z = Z_Dict.get(instance,lambda: 0)()
            Exact_Z = foo(instance)
            if solver_name in ['park', 'yuan']: # no ibound
                template_solver_filled = MEM_JOB.format(
                    use_slots = slot_assigned,
                    use_clusters = cluster_assigned,
                    init_dir=haplo1_root + '/' + solver_name + '/' + domain + '/' + instance,
                    uai_file = uai_file,
                    evid_file = evid_file,
                    vo_file = vo_file,
                    option_strings_name=''.join(command_per_solver[solver_name].split()),
                    option_strings_cmd = command_per_solver[solver_name],
                    program = '/build/' + program_name,
                    time_limit_shell = 3600,
                    time_limit = time_limit,
                    memory_limit_kb = memory_limit*1000,
                    i_bd= 0
                    )
                condor_file.write(template_solver_filled)

            else:
                for i_bd in ibounds:
                    template_solver_filled = MEM_JOB.format(
                        use_slots = slot_assigned,
                        use_clusters = cluster_assigned,
                      #  init_dir= haplo1_root + '/'  + domain + '/' + instance,
                        init_dir= '.',
                        uai_file = uai_file,
                        uai_loc = '/home/hnguyen/Experiments/'  + domain + '/' + instance + '/'+ uai_file,
                        evid_file = evid_file,
                        evid_loc = '/home/hnguyen/Experiments/'  + domain + '/' + instance + '/' + evid_file,
                        vo_file = vo_file,
                        vo_loc = '/home/hnguyen/Experiments/'  + domain + '/' + instance + '/' + vo_file,
                        option_strings_name=''.join(command_per_solver[solver_name].split()),
                        option_strings_cmd = command_per_solver[solver_name],
                        program = '/home/hnguyen/Experiments/build/' + program_name,
                        time_limit_shell = 3600,
                        time_limit = time_limit,
                        memory_limit_kb = memory_limit*1000,
                        i_bd= i_bd,
                        Z_value = Exact_Z
                        )
                    condor_file.write(template_solver_filled)

####################################################################################################################

run(solver_name='EXP_AO_8_N', time_limit=3600)