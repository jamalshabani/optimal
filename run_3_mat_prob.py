import subprocess

# Ratio Er/Es = 10.0
program_list = ["rm -rf output1",
                "python3 3_mat_prob.py -tao_monitor -tao_max_it 300 -lr 0.3 -ls 0.3 -m 'main.msh' -o 'output1' -er 1.0e-1 -es 1.0 -k 2.2e-3 -e 4.0e-3 -p 2.0",]
                #"python3 blocking_load_2_mat.py -tao_monitor -tao_max_it 1000 -l 0.10 -m '1_to_1_mesh.msh' -o '1_to_1_mesh' -er 1.0 -es 1.0e-1 -k 1.0e-5 -e 4.0e-3 -p 1.0",]
                #"python3 blocking_load_2_mat.py -tao_monitor -tao_max_it 1000 -l 0.016 -m '1_to_6_mesh.msh' -o '1_to_6_mesh' -er 1.0 -es 1.0e-1 -k 1.0e-5 -e 4.0e-3 -p 1.0"]

# Ratio Er/Es = 1.0
# program_list = ["python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l 0.25 -m '1_to_1_mesh.msh' -o '1_to_1_mesh' -er 1.0 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0",
#                 "python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l 1.0 -m '1_to_3_mesh.msh' -o '1_to_3_mesh' -er 1.0 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0",
#                 "python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l 2.0 -m '1_to_6_mesh.msh' -o '1_to_6_mesh' -er 1.0 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0"]

# Ratio Er/Es = 0.1
# program_list = ["python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l -0.005 -m '1_to_1_mesh.msh' -o '1_to_1_mesh' -er 1.0e-1 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0",
#                 "python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l -0.0001 -m '1_to_3_mesh.msh' -o '1_to_3_mesh' -er 1.0e-1 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0",
#                 "python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l -0.01 -m '1_to_6_mesh.msh' -o '1_to_6_mesh' -er 1.0e-1 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0"]

i = 1
for program in program_list:
    print("------------------------------------------------------------------------------")
    print("")
    print("Running test #{}".format(i))
    print("")
    print(program)
    print("")
    subprocess.run(program, shell = True)
    i = i + 1
