import hashlib
import json
import os
import io
from typing import Any, Dict
from Bio import Phylo


class FileCaching:
    CACHE_FILE = "./results/cache.json"
    _cache: Dict[str, str] = {}
    _is_loaded = False

    @classmethod
    def _ensure_cache_loaded(cls) -> None:
        """Ensure the cache is loaded from disk if it hasn't been already."""
        if not cls._is_loaded:
            if os.path.exists(cls.CACHE_FILE):
                try:
                    with open(cls.CACHE_FILE, "r") as f:
                        cls._cache = json.load(f)
                except json.JSONDecodeError:
                    cls._cache = {}
            cls._is_loaded = True

    @classmethod
    def _save_cache(cls) -> None:
        """Save the current cache to the JSON file."""
        with open(cls.CACHE_FILE, "w") as f:
            json.dump(cls._cache, f, indent=2)

    @staticmethod
    def _calculate_hash(file_path: str) -> str:
        """Calculate SHA-256 hash of the file content."""
        sha256_hash = hashlib.sha256()

        with open(file_path, "rb") as f:
            # Read and update hash in chunks for memory efficiency
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)

        return sha256_hash.hexdigest()

    @classmethod
    def save_genetic_tree_result(cls, input_file_path: str, output_file_path: str, aligments: Dict[str, Any], genetic_trees: Dict[str, Any]) -> None:
        """Save the genetic tree result to the cache."""
        with open(output_file_path, "w") as f:
            json.dump({"msa": aligments, "geneticTrees": genetic_trees}, f)
        cls.cache_result(input_file_path, output_file_path)

    @classmethod
    def get_cached_result_file(cls, file_path: str) -> Dict[str, Any] | None:
        """Get cached analysis result if it exists."""
        cls._ensure_cache_loaded()
        file_hash = cls._calculate_hash(file_path)
        cached_file = cls._cache.get(file_hash)
        if cached_file is None:
            return None
        result = None
        try:
            with open(cached_file, "r") as f:
                result = json.load(f)
            trees = {}
            for key, value in result["geneticTrees"].items():
                with io.StringIO(value.strip()) as file:
                    trees[key] = Phylo.read(file, format="newick")
            result["geneticTrees"] = trees
                
        except OSError:
            print("Could not get the file in the cache")
            cls._cache.pop(file_hash)
            cls._save_cache()
            return None

        return result

    @classmethod
    def cache_result(cls, file_path: str, result_file: str) -> None:
        """Cache the analysis result for the given file."""
        cls._ensure_cache_loaded()
        file_hash = cls._calculate_hash(file_path)
        cls._cache[file_hash] = result_file
        cls._save_cache()

    @classmethod
    def clear_cache(cls) -> None:
        """Clear both in-memory and file cache."""
        cls._cache = {}
        cls._is_loaded = False
        if os.path.exists(cls.CACHE_FILE):
            os.remove(cls.CACHE_FILE)
