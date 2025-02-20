class Node:
    def __init__(self, key):
        self.key = key
        self.left = None
        self.right = None


class BinaryTree:
    def __init__(self):
        # Private member variable for the root node
        self._root_node = None

    def insert(self, key):
        """Public method to insert a key into the binary tree."""
        self._root_node = self._insert_recursive(self._root_node, key)

    def _insert_recursive(self, node, key):
        """Helper method to recursively insert a key into the tree."""
        if node is None:
            return Node(key)
        if key < node.key:
            node.left = self._insert_recursive(node.left, key)
        else:
            node.right = self._insert_recursive(node.right, key)
        return node

    def inorder_traversal(self):
        """Public method to perform in-order traversal and print the tree's keys."""
        self._inorder_recursive(self._root_node)
        print()  # Newline after traversal

    def _inorder_recursive(self, node):
        """Helper method for recursive in-order traversal."""
        if node:
            self._inorder_recursive(node.left)
            print(node.key, end=" ")
            self._inorder_recursive(node.right)