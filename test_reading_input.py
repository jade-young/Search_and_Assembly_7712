"""Description: Unit tests for the reading input functions
"""
import unittest
import reading_input


class TestInput(unittest.TestCase):
    # test that the read_reads function actually returns a dictionary
    def test_read_reads(self):
        self.assertIsNotNone(reading_input.read_reads("test_fasta.txt"))

        # test that the dictionary returned by the read_reads function has the same 
        # number of keys as lines in the given file
        with open("test_fasta.txt") as file:
            file_length = len(file.readlines())
            self.assertEqual(len(reading_input.read_reads("test_fasta.txt")), file_length/2)

    def test_read_query(self):
        self.assertIsNotNone(reading_input.read_query("test_query.txt"))
        # test that the dictionary returned by the read_query function has the same
        # number of keys as lines in the given file
        with open("test_query.txt") as file:
            file_length = len(file.readlines())
            self.assertEqual(len(reading_input.read_query("test_query.txt")), file_length/2)

if __name__ == '__main__':
    unittest.main()
