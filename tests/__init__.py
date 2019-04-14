import unittest
import os
import tempfile
import shutil


class OSmiPyTestCase(unittest.TestCase):

    TEST_FILES = 'tests_files'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if kwargs.get('use_test_dir', False):
            self.tests_files_directory = os.path.join(os.path.dirname(__file__), self.TEST_FILES)
            self.assertTrue(
                os.path.exists(self.tests_files_directory),
                msg='test directory {} does not exist'.format(self.tests_files_directory))

        self.temporary_directory = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temporary_directory)

    def copy_to_temporary_directory(self, path, new_name=''):
        """Copy the content of a file to the temporary directory

        :param path: path to the file to copy
        :type path: str
        :param new_name: the new name of the file in the temporary directory (if blank, the one from path is used)
        :type new_name: str
        :rtype: str
        """

        path_in_test = os.path.join(self.tests_files_directory, path)

        if not os.path.exists(path_in_test):
            raise FileNotFoundError(path_in_test)

        if not new_name:
            new_name = os.path.basename(path)

        path_in_temp = os.path.join(self.temporary_directory, new_name)

        if os.path.exists(path_in_temp):
            raise FileExistsError(path_in_temp)

        with open(path_in_temp, 'wb') as f:
            with open(os.path.join(self.tests_files_directory, path), 'rb') as fx:
                f.write(fx.read())

        return path_in_temp
