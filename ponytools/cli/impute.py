import pkg_resources

def impute(args):
    makefile = open(
        pkg_resources.resource_filename('ponytools','data/ImputationMakefile'),mode='r'
    )
    with open('Makefile','w') as OUT:
        print(''.join(makefile.readlines()),file=OUT)


