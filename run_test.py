import subprocess

# "python3 motion_3_mat_stimulus_kI_problem.py -tao_type blmvm -tao_monitor -tao_max_it 1000 -tao_gatol 1e-7 -lr 2.0 -ls 0.14 -m 'motion_mesh.msh' -o 'motion_test1_kI' -er 1.0 -es 1.0e-2 -k 1.7e-3 -e 4.0e-3 -p 2.0 -s 0.0",
# "python3 motion_3_mat_stimulus_kI_problem.py -tao_type blmvm -tao_monitor -tao_max_it 1000 -tao_gatol 1e-7 -lr 2.0 -ls 0.13 -m 'motion_mesh.msh' -o 'motion_test2_kI' -er 1.0 -es 1.0e-2 -k 1.7e-3 -e 4.0e-3 -p 2.0 -s 0.0"

# Ratio Er/Es = 10.0
program_list = ["rm -rf test",
                "python3 test.py -tao_type pdipm -tao_monitor -tao_max_it 1000 -tao_gatol 1e-7 -lr -0.01 -ls -0.00 -m 'motion_mesh1.msh' -o 'test' -er 1.0 -es 1.0e-2 -k 1.0e-3 -e 4.0e-3 -p 2.0 -s 1.0"]


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
