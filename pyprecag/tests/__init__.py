import os
import shutil


def setup_folder(base_folder, new_folder=''):
    """Usage:
        Folder from class name:  self.TmpDir = setup_folder(__class__.__name__)
        Folder from testname:    out_dir = setup_folder(self.TmpDir, new_folder=self._testMethodName)
    """

    if new_folder in base_folder or new_folder == '':
        out_path = base_folder
    else:
        out_path = os.path.join(base_folder, new_folder)

    if os.path.exists(out_path):
        print('Folder Exists.. Deleting {}'.format(out_path))
        shutil.rmtree(out_path)
    os.makedirs(out_path)

    return out_path