from PyQt6.QtCore import QObject, pyqtSignal
from aphylogeo import utils
from aphylogeo.alignement import AlignSequences
from aphylogeo.genetic_trees import GeneticTrees
from aphylogeo.params import Params


class Worker(QObject):
    progress = pyqtSignal(int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, filepath):
        super().__init__()
        self.filepath = filepath
        self.running = True

    def run(self):
        try:
            # Step 1: Load sequences
            self.progress.emit(0)
            if not self.running:
                return
            sequenceFile = utils.loadSequenceFile(self.filepath)  # noqa: N806

            # Step 2: Align sequences
            self.progress.emit(1)
            if not self.running:
                return
            align_sequence = AlignSequences(sequenceFile)
            alignments = align_sequence.align()

            # Step 3: Generate genetic trees
            self.progress.emit(2)
            if not self.running:
                return
            geneticTrees = utils.geneticPipeline(alignments.msa)  # noqa: N806

            trees = GeneticTrees(trees_dict=geneticTrees, format="newick")

            # Step 4: Preparing results
            msa = alignments.to_dict().get("msa")

            # Step 5: Save results
            alignments.save_to_json(f"./scripts/results/aligned_{Params.reference_gene_file}.json")
            trees.save_trees_to_json("./scripts/results/geneticTrees.json")

            # Emit finished signal with the genetic trees dictionary
            result = {"msa": msa, "geneticTrees": geneticTrees}
            self.progress.emit(3)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(Exception(f"{e} consider to change the Tree Type in the alignment settings")))

    def stop(self):
        self.running = False