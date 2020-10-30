from yaml   import safe_load as yaml_safe_load
from yaml   import YAMLError
from os     import path      as os_path
from shutil import copyfile

# input files
_env_file     = 'ci/pytest/_environment.yml'
env_file      = 'ci/pytest/environment.yml'
channels_file = 'recipe/conda_channels.txt'
recipe_file   = 'recipe/meta.yaml'
bld_cfg_file  = 'recipe/conda_build_config.yaml'

# output files
cmd_file    = 'ci/pytest/cmd.sh'

def parse_meta(filename):

    recipe = ''
    with open(filename, 'r') as f:
        line = f.readline()
        while line:
            # filter all non-YAML elements
            if not line.startswith('{%') and '{{' not in line:
                recipe += line
            line = f.readline()

    requirements = []
    try:
        try: requirements += yaml_safe_load(recipe)['requirements']['host']
        except TypeError: pass
        try: requirements += yaml_safe_load(recipe)['requirements']['run']
        except TypeError: pass
        try: requirements += yaml_safe_load(recipe)['test']['requires']
        except TypeError: pass
        tests_cmd = yaml_safe_load(recipe)['test']['commands']
    except YAMLError as exc:
        print(exc)

    return requirements, tests_cmd

def parse_build_config(filename):
    with open(filename, 'r') as f:
        f_configs = f.read()
        try: configs = yaml_safe_load(f_configs)
        except TypeError: pass
    return configs


def write_dependencies(filename, requirements):

    channels = open(channels_file, 'r').read().split()

    with open(filename, 'w') as f:
        f.write('channels:\n')
        for c in channels:
            f.write('  - '+c+'\n')
        f.write('dependencies:\n')
        for req in requirements:
            f.write('  - '+req+'\n')

def write_commands(filename, cmd):
    with open(filename, 'w') as f:
        f.write(' ; '.join(tests_cmd))


requirements, tests_cmd = parse_meta(recipe_file)
write_dependencies(env_file, requirements)
write_commands(cmd_file, tests_cmd)

bld_cfgs = parse_build_config(bld_cfg_file)
for pkg in bld_cfgs:
    for ver in bld_cfgs[pkg]:
        env_file_cfg = env_file+'.'+pkg+str(ver)
        dep = pkg+'='+str(ver)
        print(dep)
        # write_dependencies(env_file_cfg, requirements+[dep])
