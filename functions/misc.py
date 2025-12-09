import sys

def package_check():
    module_names = ['astropy', 'matplotlib', 'math', 'numpy', 'pandas'];
    # module_name = 'math';
    for i in range(0, len(module_names)):
        if module_names[i] not in sys.modules:
    # txt = '{}'
            print('{} is NOT imported'.format(module_names[i]));
        else:
            print('{} is imported'.format(module_names[i]));