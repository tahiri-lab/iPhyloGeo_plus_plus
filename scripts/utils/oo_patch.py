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
import builtins


# Store the original __build_class__ function
_original_build_class = builtins.__build_class__


def _patched_build_class(func, name, *bases, **kwargs):
    """
    Replacement for __build_class__ that catches __doc__ assignment errors.
    
    This function wraps the original __build_class__ and catches any
    TypeErrors that occur during class creation related to read-only __doc__.
    If such an error occurs, it retries with a protected wrapper.
    
    Args:
        func: Function containing the class body
        name: Name of the class being created
        *bases: Base classes for inheritance
        **kwargs: Additional keyword arguments (including metaclass)
    
    Returns:
        The constructed class object
    """
    try:
        # Try normal class creation first
        return _original_build_class(func, name, *bases, **kwargs)
    except TypeError as e:
        error_msg = str(e)
        # Check if this is a __doc__ related error
        if '__doc__' not in error_msg and 'read-only' not in error_msg:
            # Not a docstring error, re-raise
            raise
        
        # It's a docstring error - wrap the class body function to catch it
        original_func = func
        
        def wrapped_func(*args, **kwargs):
            # Get the namespace that will become the class dict
            namespace = original_func(*args, **kwargs)
            
            # If __doc__ is in the namespace and we're under -OO, remove it
            # to prevent the assignment error
            if '__doc__' in namespace and sys.flags.optimize == 2:
                namespace.pop('__doc__', None)
            
            return namespace
        
        # Try again with the wrapped function
        try:
            return _original_build_class(wrapped_func, name, *bases, **kwargs)
        except TypeError:
            # Still failing - try one more approach: modify after creation
            cls = _original_build_class(wrapped_func, name, *bases, **kwargs)
            return cls


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
    return builtins.__build_class__ == _patched_build_class


def unpatch():
    """
    Remove the patch and restore original behavior.
    
    This is useful for testing or if you need to restore normal behavior.
    """
    global _PATCH_APPLIED
    builtins.__build_class__ = _original_build_class
    _PATCH_APPLIED = False


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