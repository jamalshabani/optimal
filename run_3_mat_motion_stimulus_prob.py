import subprocess

# Ratio Er/Es = 10.0
program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 3_mat_motion_stimulus_prob.py -tao_type blmvm -tao_monitor -tao_max_it 3000 -tao_gatol 1e-7 -lr 2.6 -ls 0.130 -m 'motion_mesh.msh' -o 'test1_no_stimulus' -er 1.0 -es 1.0e-2 -k 1.9e-3 -e 4.0e-3 -p 2.0 -a 0.0",
                "python3 3_mat_motion_stimulus_prob.py -tao_type blmvm -tao_monitor -tao_max_it 3000 -tao_gatol 1e-7 -lr 2.6 -ls 0.130 -m 'motion_mesh.msh' -o 'test2_no_stimulus' -er 1.0 -es 1.0e-2 -k 2.0e-3 -e 4.0e-3 -p 2.0 -a 0.0"]
                #"python3 3_mat_motion_stimulus_prob.py -tao_type blmvm -tao_monitor -tao_max_it 3000 -tao_gatol 1e-7 -lr 2.6 -ls 0.130 -m 'motion_mesh.msh' -o 'test1' -er 1.0 -es 1.0e-2 -k 2.2e-3 -e 4.0e-3 -p 2.0 -a 1.0e-3"]
                #"python3 3_mat_motion_prob.py -tao_type blmvm -tao_monitor -tao_max_it 2000 -lr 2.8 -ls 0.160 -m 'motion_mesh.msh' -o 'test2' -er 1.0 -es 1.0e-2 -k 2.2e-3 -e 4.0e-3 -p 2.0",
                #"python3 3_mat_motion_prob.py -tao_type blmvm -tao_monitor -tao_max_it 2000 -lr 3.0 -ls 0.160 -m 'motion_mesh.msh' -o 'test3' -er 1.0 -es 1.0e-2 -k 2.2e-3 -e 4.0e-3 -p 2.0"]


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
