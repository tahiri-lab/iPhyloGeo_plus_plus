"""
Patch to make ete3, and by extension aphylogeo, compatible with Python's -OO optimization flag.

The -OO flag removes docstrings and makes __doc__ attributes read-only.
This patch intercepts class creation to prevent TypeErrors when ete3
attempts to programmatically set docstrings.

Usage:
    Import this module BEFORE importing ete3 or any module that imports ete3.
    
    Example:
        import utils.oo_patch  # Must be first
        import ete3
        # or
        import utils.oo_patch
        from aphylogeo import utils
"""


import sys


class DocstringProtectedType(type):
    """
    Metaclass that silently handles __doc__ assignment errors.
    
    Under -OO, __doc__ attributes are read-only. This metaclass
    catches and suppresses TypeErrors when code attempts to assign
    to __doc__, allowing the rest of the code to execute normally.
    """
    
    def __setattr__(cls, name, value):
        if name == '__doc__':
            try:
                super().__setattr__(name, value)
            except TypeError as e:
                # Under -OO, __doc__ is read-only. Silently ignore.
                # Store in alternate attribute if needed for debugging
                if 'read-only' in str(e) or '__doc__' in str(e):
                    # Optionally store the attempted docstring elsewhere
                    super().__setattr__('_attempted_doc', value)
                else:
                    # Re-raise if it's a different TypeError
                    raise
        else:
            super().__setattr__(name, value)


def _patched_build_class(func, name, *bases, **kwargs):
    """
    Replacement for __build_class__ that injects our protective metaclass.
    
    This function is called internally by Python whenever a class is defined.
    We intercept it to create a compatible metaclass that inherits from any
    existing metaclasses while adding our __doc__ protection.
    
    Args:
        func: Function containing the class body
        name: Name of the class being created
        *bases: Base classes for inheritance
        **kwargs: Additional keyword arguments (including metaclass)
    
    Returns:
        The constructed class object
    """
    # Determine the metaclass to use
    explicit_metaclass = kwargs.get('metaclass')
    
    # Collect metaclasses from base classes
    base_metaclasses = [type(base) for base in bases if isinstance(base, type)]
    
    # Determine the most derived metaclass
    if explicit_metaclass:
        base_metas = base_metaclasses
    else:
        base_metas = base_metaclasses if base_metaclasses else [type]
        explicit_metaclass = base_metas[0] if base_metas else type
    
    # If we need protection and the metaclass isn't already protected
    if explicit_metaclass != DocstringProtectedType and not issubclass(explicit_metaclass, DocstringProtectedType):
        # Create a new metaclass that combines the existing one with our protection
        protected_metaclass = type(
            f'Protected{explicit_metaclass.__name__}',
            (DocstringProtectedType, explicit_metaclass),
            {}
        )
        kwargs['metaclass'] = protected_metaclass
    
    return _original_build_class(func, name, *bases, **kwargs)


# Store the original __build_class__ function
# __builtins__ can be a module or dict depending on context
import builtins
_original_build_class = builtins.__build_class__

# Replace it with our patched version
builtins.__build_class__ = _patched_build_class


# Flag to indicate the patch has been applied
_PATCH_APPLIED = True


def is_patch_active():
    """
    Check if the patch is currently active.
    
    Returns:
        bool: True if the patch has been applied
    """
    import builtins
    return builtins.__build_class__ == _patched_build_class


def get_patch_info():
    """
    Get information about the patch status.
    
    Returns:
        dict: Dictionary containing patch status information
    """
    return {
        'patch_applied': _PATCH_APPLIED,
        'patch_active': is_patch_active(),
        'python_optimization_level': sys.flags.optimize,
        'docstrings_stripped': sys.flags.optimize == 2,
    }