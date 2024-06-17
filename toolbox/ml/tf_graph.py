# SPDX-License-Identifier: LGPL-3.0-or-later
import os

from deepmd.env import tf


class Graph:
    """
    from toolbox.ml import Graph

    Graph("graph.pb").run()
    # then you can access the graph with localhost:6006
    """

    def __init__(self, dp_model="graph.pb") -> None:
        self.model = dp_model

    def run(self, port=6006):
        graph = self._load_graph(self.model)
        with tf.Session(graph=graph) as sess:
            writer = tf.summary.FileWriter("dp_logs", sess.graph)
            writer.close()
        os.system("tensorboard --logdir=dp_logs --port=%d" % port)

    @staticmethod
    def _load_graph(
        frozen_graph_filename, prefix: str = "load", default_tf_graph: bool = False
    ):
        """
        deepmd.infer.DeepEval._load_graph
        """
        # We load the protobuf file from the disk and parse it to retrieve the
        # unserialized graph_def
        with tf.gfile.GFile(str(frozen_graph_filename), "rb") as f:
            graph_def = tf.GraphDef()
            graph_def.ParseFromString(f.read())

            if default_tf_graph:
                tf.import_graph_def(
                    graph_def,
                    input_map=None,
                    return_elements=None,
                    name=prefix,
                    producer_op_list=None,
                )
                graph = tf.get_default_graph()
            else:
                # Then, we can use again a convenient built-in function to import
                # a graph_def into the  current default Graph
                with tf.Graph().as_default() as graph:
                    tf.import_graph_def(
                        graph_def,
                        input_map=None,
                        return_elements=None,
                        name=prefix,
                        producer_op_list=None,
                    )

            return graph
