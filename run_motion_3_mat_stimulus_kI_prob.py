import subprocess


# Ratio Er/Es = 10.0
program_list = ["rm -rf motion_test1_kI", "rm -rf motion_test2_kI",
                "python3 motion_3_mat_stimulus_kI_prob.py -tao_type blmvm -tao_monitor -tao_max_it 5000 -tao_gatol 1e-7 -lr 1.0 -ls 0.1 -m 'motion_mesh.msh' -o 'motion_test1_kI' -er 1.0 -es 1.0e-2 -k 1.9e-3 -e 4.0e-3 -p 2.0 -a 1.0 -s 1.0",
                "python3 motion_3_mat_stimulus_kI_prob.py -tao_type blmvm -tao_monitor -tao_max_it 5000 -tao_gatol 1e-7 -lr 1.0 -ls 0.1 -m 'motion_mesh.msh' -o 'motion_test2_kI' -er 1.0 -es 1.0e-2 -k 2.0e-3 -e 4.0e-3 -p 2.0 -a 1.0e-2 -s 1.0"]


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
