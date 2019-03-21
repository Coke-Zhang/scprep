import importlib

# Key:
# { module : [{submodule: [method1, class1]}, method2, class2]}
# module loads method2 and class2 when either of them is accessed
# module.submodule exists from initialization, but only loads
# module.submodule.method1 and module.submodule.method2 when
# either of them is accessed.
importspec = {
    'matplotlib': [{'colors': ['is_color_like', 'to_rgba', 'Normalize',
                               'ListedColormap', 'LinearSegmentedColormap',
                               'LogNorm', 'SymLogNorm', 'PowerNorm'],
                    'pyplot': ['scatter', 'rcParams', 'subplots',
                               'hist', 'setp', 'close', 'show'],
                    'animation': ['FuncAnimation'],
                    'cm': ['tab10'],
                    'axes': ['Axes'],
                    'lines': ['Line2D'],
                    'ticker': ['MaxNLocator'],
                    'transforms': ['ScaledTranslation']},
                   'get_backend', 'rc'],
    'mpl_toolkits': [{'mplot3d': ['Axes3D']}],
    'fcsparser': [{'api': ['ParserFeatureNotImplementedError', 'parse']}],
    'rpy2': [{'robjects': [{'numpy2ri': ['activate', 'ri2py'],
                            'packages': ['STAP'],
                            'vectors': ['ListVector'],
                            'conversion': ['ri2py']}],
              'rinterface': ['NULL', 'RRuntimeWarning']}],
    'h5py': ['File', 'Group', 'Dataset'],
    'tables': ['open_file', 'File', 'Group', 'CArray']
}


class AliasModule(object):

    def __init__(self, name, members):
        self.__module_name__ = name
        self.__module_members__ = members
        # always import these members if they exist
        self.__implicit_members__ = [
            '__version__', '__warning_registry__', '__file__',
            '__loader__', '__path__', '__doc__', '__package__']
        self.__loaded__ = False
        # create submodules
        for member in members:
            if isinstance(member, dict):
                for submodule, submembers in member.items():
                    setattr(self, submodule, AliasModule(
                        "{}.{}".format(name, submodule), submembers))

    def __getattribute__(self, attr):
        # easy access to AliasModule members to avoid recursionerror
        super_getattr = super().__getattribute__
        if attr in (super_getattr("__module_members__") +
                    super_getattr("__implicit_members__")):
            # accessing a lazy loaded member
            if not super_getattr("__loaded__"):
                # module hasn't been imported yet
                setattr(
                    self, "__module__",
                    importlib.import_module(super_getattr("__module_name__")))
            # access lazy loaded member from loaded module
            return getattr(super_getattr("__module__"), attr)
        else:
            # not a loaded member
            return super_getattr(attr)


# load required aliases into global namespace
# these can be imported as
# from scprep._lazyload import matplotlib
# plt = matplotlib.pyplot
for module, members in importspec.items():
    globals()[module] = AliasModule(module, members)
