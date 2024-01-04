# Parameters for the Musubi recheck runs.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Full list of testcases that should be in recheck:                           !
#! https://geb.inf.tu-dresden.de/collab/projects/musubi/wiki/Testcases         !
#! If you add a new testcase to recheck please update the Wiki page given above!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import os
import sys
import datetime
import shutil

from clone_and_build_function import *

templateFolder = './templates/'
machineFolder = './machines/'

date = datetime.datetime.now().strftime("%Y-%m-%d__%X")
weekday = datetime.datetime.now().strftime("%A")

# Production directory, keep the past week as history.
prod_dir = 'musubi-runs_' + weekday

run_label = 'MUSUBI'

# Cleanup production directory before using it:
shutil.rmtree(prod_dir, ignore_errors=True)
loglevel = 'INFO'

git_clone_source = 'https://github.com/apes-suite/'

from recheck import notify_list
mail_address = notify_list
smtp_server = { 'tunnel' : {'host':'robin.inf.tu-dresden.de'} }

# name of the shepherd log file
shepherd_out = 'shepherd.log'

# name of the log and rror file of the clone and build function
clone_build_out = 'clone_build.log'
clone_build_err = 'clone_build_error.log'

create_tag_on = False

try:
    # Let shepherd potentially override some settings.
    # Needed for the manual checking.
    from overrides import grep_performance, mus_rev, mus_buildvar, seeder_exe, apesFolder

except ModuleNotFoundError as error:
    apesFolder = os.path.join(os.getenv('HOME'), 'apes')
    grep_performance = True
    mus_rev = None
    mus_buildvar = 'build'

    # Use the latest seeder revision
    seeder_exe = clone_build(
            solver          = 'seeder',
            revision        = 'main',
            git_clone_source = git_clone_source+'seeder.git',
            solver_dir      = 'seeder',
            clone_build_out = clone_build_out,
            clone_build_err = clone_build_err         )

musubi_test = os.path.join(apesFolder, 'musubi', 'examples') + '/'
loris_clone_url = os.path.join(apesFolder, 'loris') + '/'

shepherd_jobs = []


musubi_exe = clone_build( solver          = 'musubi',
        revision        = mus_rev,
        variant         = mus_buildvar,
        git_clone_source = git_clone_source+'musubi.git',
        solver_dir      = 'musubi',
        clone_build_out = clone_build_out,
        clone_build_err = clone_build_err         )

musubi_omp = clone_build( solver          = 'musubi',
        git_clone_source = git_clone_source+'musubi.git',
        revision        = mus_rev,
        variant         = mus_buildvar,
        solver_dir      = 'musubi',
        confopts        = '--openmp',
        clone_build_out = clone_build_out,
        clone_build_err = clone_build_err         )

################################################################################
#                                                                              #
#                Testcases for scheme kind: fluid                              #
#                                                                              #
################################################################################
### Path to fluid testcases
fluid = musubi_test+'fluid/'


#//////////////////////////////////////////////////////////////////////////////#
#                   Testcases tree: fluid/benchmark                            #
#//////////////////////////////////////////////////////////////////////////////#
### Path to fluid benchmark testcases
fluid_bench = fluid+'benchmark/'


#------------------------------------------------------------------------------#
### start gaussianPulse testcase
testcase_path = fluid_bench+'gaussianPulse/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'gaussianPulse',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 2',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    create_dir = True,
    label = 'gaussianPulse_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/gaussianPulse_pressAlongLength_p00000_t10.001E+00.res',
    val_output_filename = 'tracking/gaussianPulse_pressAlongLength_p00000_t10.001E+00.res',
    ))
### end gaussianPulse testcase
#------------------------------------------------------------------------------#


#==============================================================================#
#                Testcases tree: fluid/benchmark/Channel2D                     #
#==============================================================================#
### Path to fluid benchmark Channel2D testcases
fluid_bench_C2D = fluid_bench+'Channel2D/'


#==============================================================================#
#        Testcases tree: fluid/benchmark/Channel2D/C2D_Simple/                 #
#==============================================================================#
### Path to fluid benchmark C2D_Simple testcases
fluid_bench_C2D_Simple = fluid_bench_C2D+'C2D_Simple/'


#------------------------------------------------------------------------------#
### start C2D_Simple_BGK
testcase_path = fluid_bench_C2D_Simple+'C2D_Simple_BGK/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_Simple_BGK',
    label = 'C2D_Simple_BGK_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 2',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Simple_BGK_seeder'],
    create_dir = False,
    label = 'C2D_Simple_BGK_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [0,4],
    val_ref_path = testcase_path+'reference/channel_pressAlongLength_p00000_t12.953E+00.res',
    val_output_filename = 'tracking/channel_pressAlongLength_p00000_t12.953E+00.res',
    ))
### end C2D_Simple_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_Simple_rBGK
testcase_path = fluid_bench_C2D_Simple+'C2D_Simple_rBGK/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_Simple_rBGK',
    label = 'C2D_Simple_rBGK_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 3',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Simple_rBGK_seeder'],
    create_dir = False,
    label = 'C2D_Simple_rBGK_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [0,4],
    val_ref_path = testcase_path+'reference/channel_pressAlongLength_p00000_t12.954E+00.res',
    val_output_filename = 'tracking/channel_pressAlongLength_p00000_t12.954E+00.res',
    ))
### end C2D_Simple_rBGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_Simple_rrBGK
testcase_path = fluid_bench_C2D_Simple+'C2D_Simple_rrBGK/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_Simple_rrBGK',
    label = 'C2D_Simple_rrBGK_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Simple_rrBGK_seeder'],
    create_dir = False,
    label = 'C2D_Simple_rrBGK_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [0,4],
    val_ref_path = testcase_path+'reference/channel_pressAlongLength_p00000_t13.337E+00.res',
    val_output_filename = 'tracking/channel_pressAlongLength_p00000_t13.337E+00.res',
    ))
### end C2D_Simple_rrBGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_Simple_prrBGK
testcase_path = fluid_bench_C2D_Simple+'C2D_Simple_prrBGK/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_Simple_prrBGK',
    label = 'C2D_Simple_prrBGK_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 1',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Simple_prrBGK_seeder'],
    create_dir = False,
    label = 'C2D_Simple_prrBGK_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [0,4],
    val_ref_path = testcase_path+'reference/channel_pressAlongLength_p00000_t9.694E+00.res',
    val_output_filename = 'tracking/channel_pressAlongLength_p00000_t9.694E+00.res',
    ))
### end C2D_Simple_prrBGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_Simple_hrrBGK
testcase_path = fluid_bench_C2D_Simple+'C2D_Simple_hrrBGK/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_Simple_hrrBGK',
    label = 'C2D_Simple_hrrBGK_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 2',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Simple_hrrBGK_seeder'],
    create_dir = False,
    label = 'C2D_Simple_hrrBGK_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [0,4],
    val_ref_path = testcase_path+'reference/channel_pressAlongLength_p00000_t13.280E+00.res',
    val_output_filename = 'tracking/channel_pressAlongLength_p00000_t13.280E+00.res',
    ))
### end C2D_Simple_hrrBGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_Simple_drtBGK
testcase_path = fluid_bench_C2D_Simple+'C2D_Simple_drtBGK/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_Simple_drtBGK',
    label = 'C2D_Simple_drtBGK_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 2',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Simple_drtBGK_seeder'],
    create_dir = False,
    label = 'C2D_Simple_drtBGK_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [0,4],
    val_ref_path = testcase_path+'reference/channel_pressAlongLength_p00000_t13.335E+00.res',
    val_output_filename = 'tracking/channel_pressAlongLength_p00000_t13.335E+00.res',
    ))
### end C2D_Simple_drtBGK
#------------------------------------------------------------------------------#


#==============================================================================#
#        Testcases tree: fluid/benchmark/Channel2D/C2D_BoundaryConditions      #
#==============================================================================#
### Path to fluid benchmark C2D_BoundaryConditions testcases
fluid_bench_C2D_BC = fluid_bench_C2D+'C2D_BoundaryConditions/'


#------------------------------------------------------------------------------#
### start C2D_BC_MfrBB_PressExpol
testcase_path = fluid_bench_C2D_BC+'C2D_BC_MfrBB_PressExpol/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_BC_MfrBB_PressExpol',
    label = 'seeder_C2D_BC_MfrBB_PressExpol',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['seeder_C2D_BC_MfrBB_PressExpol'],
    create_dir = False,
    label = 'C2D_BC_MfrBB_PressExpol_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    #  comparing difference to analytical solution
    val_ref_path = testcase_path+'reference/channel_DiffAlongHeight_p00000_t1.961E+00.res',
    val_output_filename = 'tracking/channel_DiffAlongHeight_p00000_t1.961E+00.res',
    ))
### end C2D_BC_MfrBB_PressExpol
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_BC_MfrEq_PressEq
testcase_path = fluid_bench_C2D_BC+'C2D_BC_MfrEq_PressEq/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_BC_MfrEq_PressEq',
    label = 'seeder_C2D_BC_MfrEq_PressEq',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['seeder_C2D_BC_MfrEq_PressEq'],
    create_dir = False,
    label = 'C2D_BC_MfrEq_PressEq_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    #  comparing difference to analytical solution
    val_ref_path = testcase_path+'reference/channel_DiffAlongHeight_p00000_t3.784E+00.res',
    val_output_filename = 'tracking/channel_DiffAlongHeight_p00000_t3.784E+00.res',
    ))
### end C2D_BC_MfrEq_PressEq
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_BC_PressExpol_PressExpol
testcase_path = fluid_bench_C2D_BC+'C2D_BC_PressExpol_PressExpol/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_BC_PressExpol_PressExpol',
    label = 'seeder_C2D_BC_PressExpol_PressExpol',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['seeder_C2D_BC_PressExpol_PressExpol'],
    create_dir = False,
    label = 'C2D_BC_PressExpol_PressExpol_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
#  comparing differenc to analytical soultion
    val_ref_path = testcase_path+'reference/channel_DiffAlongHeight_p00000_t5.647E+00.res',
    val_output_filename = 'tracking/channel_DiffAlongHeight_p00000_t5.647E+00.res',
    ))
### end C2D_BC_PressExpol_PressExpol
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_BC_VelBB_PressExpol
testcase_path = fluid_bench_C2D_BC+'C2D_BC_VelBB_PressExpol/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_BC_VelBB_PressExpol',
    label = 'seeder_C2D_BC_VelBB_PressExpol',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 11',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['seeder_C2D_BC_VelBB_PressExpol'],
    create_dir = False,
    label = 'C2D_BC_VelBB_PressExpol_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
#  comparing differenc to analytical soultion
    val_ref_path = testcase_path+'reference/channel_DiffAlongHeight_p00000_t2.362E+00.res',
    val_output_filename = 'tracking/channel_DiffAlongHeight_p00000_t2.362E+00.res',
    ))
### end C2D_BC_VelBB_PressExpol
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_BC_VelBFL_PressExpol
testcase_path = fluid_bench_C2D_BC+'C2D_BC_VelBFL_PressExpol/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_BC_VelBFL_PressExpol',
    label = 'seeder_C2D_BC_VelBFL_PressExpol',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 13',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['seeder_C2D_BC_VelBFL_PressExpol'],
    create_dir = False,
    label = 'C2D_BC_VelBFL_PressExpol_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
#  comparing differenc to analytical soultion
    val_ref_path = testcase_path+'reference/channel_DiffAlongHeight_p00000_t3.666E+00.res',
    val_output_filename = 'tracking/channel_DiffAlongHeight_p00000_t3.666E+00.res',
    ))
### end C2D_BC_VelBFL_PressExpol
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_BC_VelEq_PressEq
testcase_path = fluid_bench_C2D_BC+'C2D_BC_VelEq_PressEq/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_BC_VelEq_PressEq',
    label = 'seeder_C2D_BC_VelEq_PressEq',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 3',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['seeder_C2D_BC_VelEq_PressEq'],
    create_dir = False,
    label = 'C2D_BC_VelEq_PressEq_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
#  comparing differenc to analytical soultion
    val_ref_path = testcase_path+'reference/channel_DiffAlongHeight_p00000_t2.664E+00.res',
    val_output_filename = 'tracking/channel_DiffAlongHeight_p00000_t2.664E+00.res',
    ))
### end C2D_BC_VelEq_PressEq
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C2D_BC_VelNonEqExpol_PressNonEqExpol
testcase_path = fluid_bench_C2D_BC+'C2D_BC_VelNonEqExpol_PressNonEqExpol/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'C2D_BC_VelNonEqExpol_PressNonEqExpol',
    label = 'seeder_C2D_BC_VelNonEqExpol_PressNonEqExpol',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 5',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['seeder_C2D_BC_VelNonEqExpol_PressNonEqExpol'],
    create_dir = False,
    label = 'C2D_BC_VelNonEqExpol_PressNonEqExpol_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
#  comparing differenc to analytical soultion
    val_ref_path = testcase_path+'reference/channel_DiffAlongHeight_p00000_t2.451E+00.res',
    val_output_filename = 'tracking/channel_DiffAlongHeight_p00000_t2.451E+00.res',
    ))
### end C2D_BC_VelNonEqExpol_PressNonEqExpol
#------------------------------------------------------------------------------#


#==============================================================================#
#        Testcases tree: fluid/benchmark/Channel3D/C3D_Sphere_MultiLevel_LES/  #
#==============================================================================#
### Path to fluid benchmark C3D_Sphere_MultiLevel_LES testcases
fluid_bench_C3D_Sphere = fluid_bench+'Channel3D/C3D_Sphere_MultiLevel_LES/'

#==============================================================================#
#   Testcases tree: fluid/benchmark/Channel3D/C3D_Sphere_MultiLevel_LES/D3Q19  #
#==============================================================================#
fluid_bench_C3D_Sphere_QX = fluid_bench_C3D_Sphere + 'D3Q19' + '/'
label_prefix = 'C3D_Sph_ML_LES_' + 'D19_'
    
#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D19_DRT_BGK
my_collision = 'DRT_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D19_DRT_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D19_HRR_BGK
my_collision = 'HRR_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D19_HRR_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D19_PRR_BGK
my_collision = 'PRR_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D19_PRR_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D19_R_BGK
my_collision = 'R_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D19_R_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D19_RR_BGK
my_collision = 'RR_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D19_RR_BGK
#------------------------------------------------------------------------------#

#==============================================================================#
#   Testcases tree: fluid/benchmark/Channel3D/C3D_Sphere_MultiLevel_LES/D3Q27  #
#==============================================================================#
fluid_bench_C3D_Sphere_QX = fluid_bench_C3D_Sphere + 'D3Q27' + '/'
label_prefix = 'C3D_Sph_ML_LES_' + 'D27_'

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D27_CUM17
my_collision = 'CUM17'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end CUM17
#------------------------------------------------------------------------------#
    
#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D27_DRT_BGK
my_collision = 'DRT_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D27_DRT_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D27_HRR_BGK
my_collision = 'HRR_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D27_HRR_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D27_PRR_BGK
my_collision = 'PRR_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D27_PRR_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D27_R_BGK
my_collision = 'R_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D27_R_BGK
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start C3D_Sph_ML_LES_D27_RR_BGK
my_collision = 'RR_BGK'
my_label_prefix = label_prefix + my_collision
testcase_path = fluid_bench_C3D_Sphere_QX + my_collision + '/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = my_label_prefix,
    label = my_label_prefix + '_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = [my_label_prefix + '_seeder'],
    create_dir = False,
    label = my_label_prefix + '_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel3D_sphere_force_p00000.res',
    val_output_filename = 'tracking/channel3D_sphere_force_p00000.res',
    ))
### end C3D_Sph_ML_LES_D27_RR_BGK
#------------------------------------------------------------------------------#

#==============================================================================#
#                Testcases tree: fluid/benchmark/TurbulentChannel              #
#==============================================================================#
### Path to fluid benchmark turbulent channel testcases
fluid_bench_TC = fluid_bench+'TurbulentChannel/'

#==============================================================================#
#     Testcases tree: fluid/benchmark/TurbulentChannel/TC_SingleLevel          #
#==============================================================================#
### Path to fluid benchmark turbulent channel single level testcases
fluid_bench_TC_SL = fluid_bench_TC+'TC_SingleLevel/'

#------------------------------------------------------------------------------#
### start turbulent channel single level musker fixed point
testcase_path = fluid_bench_TC_SL+'TC_SL_MuskerFixedPoint/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'TC_SL_MuskerFP',
    label = 'TC_SL_MuskerFP_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['TC_SL_MuskerFP_seeder'],
    create_dir = False,
    label = 'TC_SL_MuskerFP_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/Channel_meanVel_p00000.res',
    val_output_filename = 'tracking/Channel_meanVel_p00000.res',
    ))
### end Turbulent channel single level musker fixed point
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start turbulent channel single level musker with Newton
testcase_path = fluid_bench_TC_SL+'TC_SL_MuskerNewton/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'TC_SL_MuskerNewton',
    label = 'TC_SL_MuskerNewton_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['TC_SL_MuskerNewton_seeder'],
    create_dir = False,
    label = 'TC_SL_MusNewton_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/Channel_meanVel_p00000.res',
    val_output_filename = 'tracking/Channel_meanVel_p00000.res',
    ))
### end Turbulent channel single level musker with Newton
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start turbulent channel single level Reichardt with fixed point
testcase_path = fluid_bench_TC_SL+'TC_SL_ReichardtFixedPoint/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'TC_SL_ReichFP',
    label = 'TC_SL_ReichFP_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['TC_SL_ReichFP_seeder'],
    create_dir = False,
    label = 'TC_SL_ReichFP_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/Channel_meanVel_p00000.res',
    val_output_filename = 'tracking/Channel_meanVel_p00000.res',
    ))
### end Turbulent channel single level Reichart with fixed point
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start turbulent channel single level explicit Power-Law
testcase_path = fluid_bench_TC_SL+'TC_SL_PowerLaw/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'TC_SL_PowerLaw',
    label = 'TC_SL_PowerLaw_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['TC_SL_PowerLaw_seeder'],
    create_dir = False,
    label = 'TC_SL_PowerLaw_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/Channel_meanVel_p00000.res',
    val_output_filename = 'tracking/Channel_meanVel_p00000.res',
    ))
### end Turbulent channel single level Power-Law
#------------------------------------------------------------------------------#

### end Turbulent channel single level
#------------------------------------------------------------------------------#


#==============================================================================#
#                Testcases tree: fluid/benchmark/absorbingLayer                #
#==============================================================================#
### Path to fluid benchmark absorbingLayer testcases
fluid_bench_ABS = fluid_bench+'absorbingLayer/'

#==============================================================================#
#        Testcases tree: fluid/benchmark/absorbingLayer/acousticPulse2D/       #
#==============================================================================#
### Path to fluid benchmark acousticPulse2D testcases
fluid_bench_ABS_Pulse = fluid_bench_ABS+'acousticPulse2D/'


#------------------------------------------------------------------------------#
### start absLayerRadial
testcase_path = fluid_bench_ABS_Pulse+'absLayerRadial/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'ABS_Pulse_Radial',
    label = 'ABS_Pulse_Radial_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['ABS_Pulse_Radial_seeder'],
    create_dir = False,
    label = 'ABS_Pulse_Radial_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [1,2],
    val_ref_path = testcase_path+'reference/pulse2D_probe_p00000.res',
    val_output_filename = 'tracking/pulse2D_probe_p00000.res',
    ))
### end absLayerRadial
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start absLayerPlane
testcase_path = fluid_bench_ABS_Pulse+'absLayerPlane/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'ABS_Pulse_Plane',
    label = 'ABS_Pulse_Plane_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['ABS_Pulse_Plane_seeder'],
    create_dir = False,
    label = 'ABS_Pulse_Plane_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [1,2],
    val_ref_path = testcase_path+'reference/pulse2D_probe_p00000.res',
    val_output_filename = 'tracking/pulse2D_probe_p00000.res',
    ))
### end absLayerPlane
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start absLayerBox
testcase_path = fluid_bench_ABS_Pulse+'absLayerBox/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'ABS_Pulse_Box',
    label = 'ABS_Pulse_Box_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['ABS_Pulse_Box_seeder'],
    create_dir = False,
    label = 'ABS_Pulse_Box_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [1,2],
    val_ref_path = testcase_path+'reference/pulse2D_probe_p00000.res',
    val_output_filename = 'tracking/pulse2D_probe_p00000.res',
    ))
### end absLayerBox
#------------------------------------------------------------------------------#

#==============================================================================#
#        Testcases tree: fluid/benchmark/absorbingLayer/acousticPulse3D/       #
#==============================================================================#
### Path to fluid benchmark acousticPulse3D testcases
fluid_bench_ABS_Pulse3D = fluid_bench_ABS+'acousticPulse3D/'

#------------------------------------------------------------------------------#
### start absLayerBox3D
testcase_path = fluid_bench_ABS_Pulse3D+'absLayerBox3D/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'ABS_Pulse_Box3D',
    label = 'ABS_Pulse_Box3D_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['ABS_Pulse_Box3D_seeder'],
    create_dir = False,
    label = 'ABS_Pulse_Box3D_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [1,2],
    val_ref_path = testcase_path+'reference/pulse3D_probe_p00000.res',
    val_output_filename = 'tracking/pulse3D_probe_p00000.res',
    ))
### end absLayerRadial
#------------------------------------------------------------------------------#

#==============================================================================#
#        Testcases tree: fluid/benchmark/absorbingLayer/acousticLineSource2D/  #
#==============================================================================#
### Path to fluid benchmark absLayer acousticLineSource2D testcases
fluid_bench_ABS_LineSrc2D = fluid_bench_ABS+'acousticLineSource2D/'

#------------------------------------------------------------------------------#
### start absLayerBox3D
testcase_path = fluid_bench_ABS_LineSrc2D
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'ABS_LineSrc2D',
    label = 'ABS_LineSrc2D_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 2',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['ABS_LineSrc2D_seeder'],
    create_dir = False,
    label = 'ABS_LineSrc2D_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [3,4],
    val_ref_path = testcase_path+'reference/acousticLineSrc_line_p00001_t13.988E-03.res',
    val_output_filename = 'tracking/acousticLineSrc_line_p00001_t13.988E-03.res',
    ))
### end acousticLineSource2D
#------------------------------------------------------------------------------#


#==============================================================================#
#        Testcases tree: fluid/benchmark/absorbingLayer/acousticCylinder2D/  #
#==============================================================================#
### Path to fluid benchmark absLayer acousticCylinder2D testcases
fluid_bench_ABS_cyl2D = fluid_bench_ABS+'acousticCylinder2D/'

#------------------------------------------------------------------------------#
### start absLayerBox3D
testcase_path = fluid_bench_ABS_cyl2D
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'ABS_cyl2D',
    label = 'ABS_cyl2D_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','tracking_vtk','restart'],
    depend = ['ABS_cyl2D_seeder'],
    create_dir = False,
    label = 'ABS_cyl2D_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel_probeAt90deg_p00000.res',
    val_output_filename = 'tracking/channel_probeAt90deg_p00000.res',
    ))
### end acousticCylinder2D
#------------------------------------------------------------------------------#


#==============================================================================#
#        Testcases tree: fluid/benchmark/trackOutside  #
#==============================================================================#
### Path to fluid benchmark trackOutside
fluid_bench_trackout = fluid_bench+'trackOutside/'

#------------------------------------------------------------------------------#
### start trackOutside
testcase_path = fluid_bench_trackout
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'TO_cyl2D',
    label = 'TO_cyl2D_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','tracking_vtk','restart'],
    depend = ['TO_cyl2D_seeder'],
    create_dir = False,
    label = 'TO_cyl2D_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel_probeAt0deg_p00000.res',
    val_output_filename = 'tracking/channel_probeAt0deg_p00000.res',
    ))
### end trackOutside
#------------------------------------------------------------------------------#


#==============================================================================#
#                Testcases tree: fluid/benchmark/viscousSpongeLayer            #
#==============================================================================#
### Path to fluid benchmark viscosus sponge layer testcases
fluid_bench_VSL = fluid_bench+'viscousSpongeLayer/'

#------------------------------------------------------------------------------#
### start VSL_Radial
testcase_path = fluid_bench_VSL+'VSL_Radial/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'VSL_Radial',
    label = 'VSL_Radial_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart', 'tracking_vtk'],
    depend = ['VSL_Radial_seeder'],
    create_dir = False,
    label = 'VSL_Radial_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [1,2],
    val_ref_path = testcase_path+'reference/channel_cyl_force_p00000.res',
    val_output_filename = 'tracking/channel_cyl_force_p00000.res',
    ))
### end VSL_Radial
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start VSL_Box2D
testcase_path = fluid_bench_VSL+'VSL_Box2D/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'VSL_Box2D',
    label = 'VSL_Box2D_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['VSL_Box2D_seeder'],
    create_dir = False,
    label = 'VSL_Box2D_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [1,2],
    val_ref_path = testcase_path+'reference/pulse2D_probe_p00000.res',
    val_output_filename = 'tracking/pulse2D_probe_p00000.res',
    ))
### end VSL_Box2D
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start VSL_Box3D
testcase_path = fluid_bench_VSL+'VSL_Box3D/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'VSL_Box3D',
    label = 'VSL_Box3D_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['VSL_Box3D_seeder'],
    create_dir = False,
    label = 'VSL_Box3D_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [1,2],
    val_ref_path = testcase_path+'reference/pulse3D_probe_p00000.res',
    val_output_filename = 'tracking/pulse3D_probe_p00000.res',
    ))
### end VSL_Box3D
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
### start VSL_Plane
testcase_path = fluid_bench_VSL+'VSL_Plane/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'VSL_Plane',
    label = 'VSL_Plane_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['VSL_Plane_seeder'],
    create_dir = False,
    label = 'VSL_Plane_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [3,4],
    val_ref_path = testcase_path+'reference/acousticLineSrc_line_p00003_t13.988E-03.res',
    val_output_filename = 'tracking/acousticLineSrc_line_p00003_t13.988E-03.res',
    ))
### end VSL_Plane
#------------------------------------------------------------------------------#


################################################################################
#                                                                              #
#                Testcases for scheme kind: fluid_incompressible               #
#                                                                              #
################################################################################
### Path to fluid incompressible testcases
fluid_incomp = musubi_test+'fluid_incompressible/'


#//////////////////////////////////////////////////////////////////////////////#
#               Testcases tree: fluid_incompressible/benchmark                 #
#//////////////////////////////////////////////////////////////////////////////#
### Path to fluid_incompressible benchmark testcases
fluid_incomp_bench = fluid_incomp+'benchmark/'


#------------------------------------------------------------------------------#
### start gaussianPulse_incomp
testcase_path = fluid_incomp_bench+'gaussianPulse/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'gaussianPulse_incomp',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 2',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    create_dir = True,
    label = 'gaussianPulse_incomp',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/gaussianPulse_pressAlongLength_p00000_t10.001E+00.res',
    val_output_filename = 'tracking/gaussianPulse_pressAlongLength_p00000_t10.001E+00.res',
    ))
### end gaussianPulse_incomp
#------------------------------------------------------------------------------#


#==============================================================================#
#        Testcases tree: fluid_incompressible/benchmark/TaylorGreenVortex      #
#==============================================================================#
### Path to fluid_incompressible benchmark TaylorGreenVortex testcases
fluid_incomp_bench_TGV = fluid_incomp_bench +'TaylorGreenVortex/'


#==============================================================================#
#  Testcases tree: fluid_incompressible/benchmark/TaylorGreenVortex/TGV_Simple #
#==============================================================================#
### Path to fluid_incompressible benchmark TaylorGreenVortex testcases
fluid_incomp_bench_TGV_Simple = fluid_incomp_bench_TGV +'TGV_Simple/'


#------------------------------------------------------------------------------#
### start TGV_Simple_Re800
testcase_path = fluid_incomp_bench_TGV_Simple+'TGV_Simple_Re800/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'TGV_Simple_Re800',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    create_dir = True,
    label = 'TGV_Simple_Re800',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/TGV_Simple_Re800_probeAtCenter_p00000.res',
    val_output_filename = 'tracking/TGV_Simple_Re800_probeAtCenter_p00000.res',
    ))
### end TGV_Simple_Re800
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start TGV_Simple_Re1600
testcase_path = fluid_incomp_bench_TGV_Simple+'TGV_Simple_Re1600/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'TGV_Simple_Re1600',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart','vtkfiles'],
    create_dir = True,
    label = 'TGV_Simple_Re1600',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/TGV_Simple_Re1600_kE_all_p00000.res',
    val_output_filename = 'tracking/TGV_Simple_Re1600_kE_all_p00000.res',
    ))
### end TGV_SimpleRe_1600
#------------------------------------------------------------------------------#

#==============================================================================#
#  Testcases tree: fluid_incompressible/benchmark/TaylorGreenVortex/TGV_LES #
#==============================================================================#
### Path to fluid_incompressible benchmark TaylorGreenVortex testcases
fluid_incomp_bench_TGV_LES = fluid_incomp_bench_TGV +'TGV_LES/'


#------------------------------------------------------------------------------#
### start TGV_LES_WALE
testcase_path = fluid_incomp_bench_TGV_LES+'TGV_WALE/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'TGV_LES_WALE',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart','vtkfiles'],
    create_dir = True,
    label = 'TGV_LES_WALE',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/TGV_LES_WALE_kE_all_p00000.res',
    val_output_filename = 'tracking/TGV_LES_WALE_kE_all_p00000.res',
    ))
### end TGV_LES_WALE
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start TGV_LES_Vreman
testcase_path = fluid_incomp_bench_TGV_LES+'TGV_Vreman/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'TGV_LES_Vreman',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart','vtkfiles'],
    create_dir = True,
    label = 'TGV_LES_Vreman',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/TGV_LES_Vreman_kE_all_p00000.res',
    val_output_filename = 'tracking/TGV_LES_Vreman_kE_all_p00000.res',
    ))
### end TGV_LES_Vreman
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start TGV_LES_SmagPDF
testcase_path = fluid_incomp_bench_TGV_LES+'TGV_SmagPDF/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'TGV_LES_SmagPDF',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart','vtkfiles'],
    create_dir = True,
    label = 'TGV_LES_SmagPDF',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/TGV_LES_SmagPDF_kE_all_p00000.res',
    val_output_filename = 'tracking/TGV_LES_SmagPDF_kE_all_p00000.res',
    ))
### end TGV_LES_SmagPDF
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start TGV_LES_SmagGradU
testcase_path = fluid_incomp_bench_TGV_LES+'TGV_SmagGradU/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    prefix = 'TGV_LES_SmagGradU',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart','vtkfiles'],
    create_dir = True,
    label = 'TGV_LES_SmagGradU',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/TGV_LES_SmagGradU_kE_all_p00000.res',
    val_output_filename = 'tracking/TGV_LES_SmagGradU_kE_all_p00000.res',
    ))
### end TGV_LES_SmagGradU
#------------------------------------------------------------------------------#


#==============================================================================#
#            Testcases tree: fluid_incompressible/benchmark/Channel2D          #
#==============================================================================#
### Path to fluid_incompressible benchmark Channel2D testcases
fluid_incomp_bench_C2D = fluid_incomp_bench+'Channel2D/'


#------------------------------------------------------------------------------#
### start Channel 2D with cylinder single level
testcase_path = fluid_incomp_bench_C2D+'C2D_Cylinder_SingleLevel/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'C2D_Cyl_SL_Incomp',
    label = 'C2D_Cyl_SL_Incomp_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 4',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Cyl_SL_Incomp_seeder'],
    create_dir = False,
    label = 'C2D_Cyl_SL_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel_Cp_p00000_t7.500E+00.res',
    val_output_filename = 'tracking/channel_Cp_p00000_t7.500E+00.res',
    ))
### end Channel 2D with cylinder single level
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start Channel 2D with cylinder multilevel
testcase_path = fluid_incomp_bench_C2D+'C2D_Cylinder_MultiLevel/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'C2D_Cyl_ML_Incomp',
    label = 'C2D_Cyl_ML_Incomp_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C2D_Cyl_ML_Incomp_seeder'],
    create_dir = False,
    label = 'C2D_Cyl_ML_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel_cyl_force_p00000.res',
    val_output_filename = 'tracking/channel_cyl_force_p00000.res',
    ))
### end Channel 2D with cylinder multilevel
#------------------------------------------------------------------------------#


#==============================================================================#
#            Testcases tree: fluid_incompressible/benchmark/Channel3D          #
#==============================================================================#
### Path to fluid_incompressible benchmark Channel3D testcases
fluid_incomp_bench_C3D = fluid_incomp_bench+'Channel3D/'


#------------------------------------------------------------------------------#
### start Channel 3D simple with nonEqExpol BCs
testcase_path = fluid_incomp_bench_C3D+'C3D_Simple/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'C3D_Simple',
    label = 'C3D_Simple_Incomp_seeder',
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 12',
    additional_params = dict(testcase_path = testcase_path),
    create_subdir = ['tracking','restart'],
    depend = ['C3D_Simple_Incomp_seeder'],
    create_dir = False,
    label = 'C3D_Simple_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/C3D_Simple_velAlongHeight_p00000_t2.001E+00.res',
    val_output_filename = 'tracking/C3D_Simple_velAlongHeight_p00000_t2.001E+00.res',
    ))
### end Channel 3D simple with nonEqExpol BCs
#------------------------------------------------------------------------------#

#==============================================================================#
#        Testcases tree: fluid_incompressible/benchmark/ConcentricCylinders
#==============================================================================#
### Path to fluid_incompressible benchmark ConcentricCylinders testcases
fluid_incomp_bench_COC = fluid_incomp_bench+'ConcentricCylinders/'


#------------------------------------------------------------------------------#
### start Concentric cylinder coeutte flow 2D fluid_incompressible model
### for singlelevel
testcase_path = fluid_incomp_bench_COC+'COC_CoeutteFlow/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = fluid_incomp_bench_COC),
    create_subdir = ['mesh'],
    prefix = 'COC_CouetteFlow_Incomp',
    label = 'COC_CouetteFlow_Incomp_seeder',
    attachment = True,
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['COC_CouetteFlow_Incomp_seeder'],
    create_dir=False,
    label = 'COC_CouetteFlow_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [3,4,6],
    val_ref_path = testcase_path+'reference/concentricCylinder_line_p00000_t1.443E+00.res',
    val_output_filename = 'tracking/concentricCylinder_line_p00000_t1.443E+00.res',
    ))
### end ConcentricCylinder 2D fluid_incompressible model
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start Concentric cylinder coeutte flow 2D fluid_incompressible model
### for multilevel
testcase_path = fluid_incomp_bench_COC+'COC_CoeutteFlow_MultiLevel/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = fluid_incomp_bench_COC),
    create_subdir = ['mesh'],
    prefix = 'COC_CouetteFlow_ML_Incomp',
    label = 'COC_CouetteFlow_ML_Incomp_seeder',
    attachment = True,
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['COC_CouetteFlow_ML_Incomp_seeder'],
    create_dir=False,
    label = 'COC_CouetteFlow_ML_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    position = [3,4,6],
    val_ref_path = testcase_path+'reference/concentricCylinder_line_p00000_t1.488E+00.res',
    val_output_filename = 'tracking/concentricCylinder_line_p00000_t1.488E+00.res',
    ))
### end ConcentricCylinder 2D fluid_incompressible model
#------------------------------------------------------------------------------#




#==============================================================================#
#        Testcases tree: fluid_incompressible/benchmark/LidDrivenCavity        #
#==============================================================================#
### Path to fluid_incompressible benchmark LidDrivenCavity testcases
fluid_incomp_bench_LDC = fluid_incomp_bench+'LidDrivenCavity/'


#------------------------------------------------------------------------------#
### start Lid-Driven Cavity 2D fluid_incompressible model
testcase_path = fluid_incomp_bench_LDC+'LDC_Simple/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    create_subdir = ['mesh'],
    prefix = 'LDC_Simple_Incomp',
    label = 'LDC_Simple_Incomp_seeder',
    attachment = True,
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['LDC_Simple_Incomp_seeder'],
    create_dir=False,
    label = 'LDC_Simple_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'tracking/ref_lidcavity_probe_p00000.res',
    val_output_filename = 'tracking/lidcavity_probe_p00000.res',
    ))
### end Lid-Driven Cavity 2D fluid_incompressible model
#------------------------------------------------------------------------------#


#==============================================================================#
#              Testcases tree: fluid_incompressible/benchmark/Pipe             #
#==============================================================================#
### Path to fluid_incompressible benchmark Pipe testcases
fluid_incomp_bench_pipe = fluid_incomp_bench+'Pipe/'


#------------------------------------------------------------------------------#
### start Pipe flow with pressure-pressure BC fluid_incompressible model
testcase_path = fluid_incomp_bench_pipe+'PIP_Simple/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'PIP_Simple_Incomp',
    label = 'PIP_Simple_Incomp_seeder',
    attachment = True,
    ))

shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['PIP_Simple_Incomp_seeder'],
    create_dir=False,
    label = 'PIP_Simple_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/pipe_probeAtCenter_p00000.res',
    val_output_filename = 'tracking/pipe_probeAtCenter_p00000.res',
    ))
### end Pipe flow with pressure-pressure BC fluid_incompressible model
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start Pipe flow with external force fluid_incompressible model
testcase_path = fluid_incomp_bench_pipe+'PIP_Force/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = fluid_incomp_bench_pipe +'PIP_Simple/'),
    create_subdir = ['mesh'],
    prefix = 'PIP_Force_Incomp',
    label = 'PIP_Force_Incomp_seeder',
    attachment = True,
    ))

shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['PIP_Force_Incomp_seeder'],
    create_dir=False,
    label = 'PIP_Force_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/pipe_probeAtCenter_p00000.res',
    val_output_filename = 'tracking/pipe_probeAtCenter_p00000.res',
    ))
### end Pipe flow with external force fluid_incompressible model
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start Pipe flow with multilevel fluid_incompressible model
testcase_path = fluid_incomp_bench_pipe+'PIP_MultiLevel/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = fluid_incomp_bench_pipe +'PIP_Simple/'),
    create_subdir = ['mesh'],
    prefix = 'PIP_ML_Incomp',
    label = 'PIP_ML_Incomp_seeder',
    attachment = True,
    ))

shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['PIP_ML_Incomp_seeder'],
    create_dir=False,
    label = 'PIP_ML_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/pipe_probeAtCenter_p00000.res',
    val_output_filename = 'tracking/pipe_probeAtCenter_p00000.res',
    ))
### end Pipe flow with multilevel fluid_incompressible model
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start Pipe flow with multilevel and LES fluid_incompressible model
testcase_path = fluid_incomp_bench_pipe+'PIP_LES/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = fluid_incomp_bench_pipe +'PIP_Simple/'),
    create_subdir = ['mesh'],
    prefix = 'PIP_LES_Incomp',
    label = 'PIP_LES_Incomp_seeder',
    attachment = True,
    ))

shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['PIP_LES_Incomp_seeder'],
    create_dir=False,
    label = 'PIP_LES_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/pipeLES_probeAtCenter_p00000.res',
    val_output_filename = 'tracking/pipeLES_probeAtCenter_p00000.res',
    ))
### end Pipe flow with multilevel fluid_incompressible model
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
### start split pipe fluid_incompressible model
testcase_path = fluid_incomp_bench_pipe+'PIP_Split/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = fluid_incomp_bench_pipe +'PIP_Split/'),
    create_subdir = ['mesh'],
    prefix = 'PIP_Split_Incomp',
    label = 'PIP_Split_Incomp_seeder',
    attachment = True,
    ))

shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['PIP_Split_Incomp_seeder'],
    create_dir=False,
    label = 'PIP_Split_Incomp_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/pipe_probeAtCenter_p00000.res',
    val_output_filename = 'tracking/pipe_probeAtCenter_p00000.res',
    ))
### end split pipe fluid_incompressible model
#------------------------------------------------------------------------------#

################################################################################
#                                                                              #
#                       Testcases for tutorial                                 #
#                                                                              #
################################################################################
tutorial_test = musubi_test + 'tutorials/tutorial_cases/'

#------------------------------------------------------------------------------#
### start tutorial Pip_force testcase
testcase_path = tutorial_test +'tutorial_PIP_Force/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'tut_pip',
    label = 'tut_pip_seeder',
    attachment = True,
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 2',
    create_subdir = ['tracking','restart'],
    depend = ['tut_pip_seeder'],
    create_dir=False,
    label = 'tut_pip_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'ref/pipe_velAlongHeight_p00000_t1.000E+00.res',
    val_output_filename = 'tracking/pipe_velAlongHeight_p00000_t1.000E+00.res',
    ))
### end tutorial Pip_force testcase
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start tutorial gaussian_pulse testcase
testcase_path = tutorial_test +'tutorial_gaussian_pulse/'
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun -np 1',
    create_subdir = ['tracking','restart'],
    create_dir=True,
    label = 'tut_gauss_pulse',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'ref/Gausspulse_track_pressure_p00000.res',
    val_output_filename = 'tracking/Gausspulse_track_pressure_p00000.res',
    ))
### end tutorial gaussian_pulse testcase
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
### start tutorial channel testcase
testcase_path = tutorial_test +'tutorial_channelGeneric/'
shepherd_jobs.append(dict(executable = seeder_exe,
    template=testcase_path+'seeder.lua',
    extension='lua',
    run_exec = True,
    additional_params = dict(stl_path = testcase_path),
    create_subdir = ['mesh'],
    prefix = 'tut_channel',
    label = 'tut_channel_seeder',
    attachment = True,
    ))
shepherd_jobs.append(dict(executable = musubi_exe,
    solver_name = 'musubi',
    template=testcase_path+'musubi.lua',
    extension='lua',
    run_exec = True,
    run_command = 'mpirun --oversubscribe -np 8',
    create_subdir = ['tracking','restart'],
    depend = ['tut_channel_seeder'],
    create_dir=False,
    label = 'tut_channel_musubi',
    attachment = True,
    validation = True,
    val_method = 'difference',
    val_ref_path = testcase_path+'reference/channel_probeAtCenter_p00000.res',
    val_output_filename = 'tracking/channel_probeAtCenter_p00000.res',
    ))
### end tutorial channel testcase
#------------------------------------------------------------------------------#
