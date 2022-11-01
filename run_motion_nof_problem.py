import subprocess


# Ratio Er/Es = 10.0
program_list = ["rm -rf motion_nof1", "rm -rf motion_nof2", "rm -rf motion_nof3",
                "python3 motion_nof_problem.py -tao_type blmvm -tao_monitor -tao_max_it 200 -tao_gatol 1e-7 -lr 0.35 -ls 0.35 -m 'motion_nof_mesh1.msh' -o 'motion_nof1' -er 1.0 -es 1.0e-2 -k 1.7e-3 -e 0.5e-3 -p 2.0 -s 0.5"]


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
