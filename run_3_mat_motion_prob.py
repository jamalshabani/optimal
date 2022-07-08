import subprocess

# Ratio Er/Es = 10.0
program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 3_mat_motion_prob.py -tao_type blmvm -tao_monitor -tao_max_it 2000 -lr 2.8 -ls 0.160 -m 'motion_mesh.msh' -o 'test1' -er 1.0 -es 1.0e-2 -k 2.2e-3 -e 4.0e-3 -p 2.0",
                "python3 3_mat_motion_prob.py -tao_type blmvm -tao_monitor -tao_max_it 2000 -lr 3.0 -ls 0.160 -m 'motion_mesh.msh' -o 'test2' -er 1.0 -es 1.0e-2 -k 2.2e-3 -e 4.0e-3 -p 2.0",
                "python3 3_mat_motion_prob.py -tao_type blmvm -tao_monitor -tao_max_it 2000 -lr 3.2 -ls 0.160 -m 'motion_mesh.msh' -o 'test3' -er 1.0 -es 1.0e-2 -k 2.2e-3 -e 4.0e-3 -p 2.0"]

# Ratio Er/Es = 1.0
# program_list = ["python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l 0.25 -m '1_to_1_mesh.msh' -o '1_to_1_mesh' -er 1.0 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0",
#                 "python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l 1.0 -m '1_to_3_mesh.msh' -o '1_to_3_mesh' -er 1.0 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0",
#                 "python3 2_mat_res_prob.py -tao_monitor -tao_max_it 200 -l 2.0 -m '1_to_6_mesh.msh' -o '1_to_6_mesh' -er 1.0 -es 1.0 -k 1.0e-5 -e 4.0e-3 -p 1.0"]

# Ratio Er/Es = 0.
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
