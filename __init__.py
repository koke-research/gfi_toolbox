# __init__.py
def classFactory(iface):
    # Import our main plugin class
    from .gfi_main import GFIPlugin
    return GFIPlugin(iface)