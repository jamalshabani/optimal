import subprocess

# Ratio Er/Es = 10.0
program_list = ["rm -rf test1", "rm -rf test2", "rm -rf test3",
                "python3 3_mat_motion_prob.py -tao_type blmvm -tao_monitor -tao_max_it 2000 -tao_gatol 1e-7 -tao_grtol 1e-7 -tao_gttol 1e-7 -lr 2.6 -ls 0.160 -m 'motion_mesh2.msh' -o 'test1' -er 1.0 -es 1.0e-2 -k 2.2e-3 -e 4.0e-3 -p 2.0"]
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
