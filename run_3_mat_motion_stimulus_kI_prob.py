import subprocess

# Ratio Er/Es = 10.0
program_list = ["python3 3_mat_motion_stimulus_kI_prob.py -tao_type blmvm -tao_monitor -tao_max_it 5000 -tao_gatol 1e-7 -lr 2.6 -ls 0.130 -m 'motion_mesh.msh' -o 'test1_kI' -er 1.0 -es 1.0e-2 -k 1.9e-3 -e 4.0e-3 -p 2.0 -a 1.0 -s 1.0",
                "python3 3_mat_motion_stimulus_kI_prob.py -tao_type blmvm -tao_monitor -tao_max_it 5000 -tao_gatol 1e-7 -lr 2.6 -ls 0.130 -m 'motion_mesh.msh' -o 'test2_kI' -er 1.0 -es 1.0e-2 -k 2.0e-3 -e 4.0e-3 -p 2.0 -a 1.0e-2 -s 1.0"]


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
