from PyQt5 import uic, QtWidgets
import sys

import os

path = os.path.dirname(os.path.abspath('__file__'))


class Shell(*uic.loadUiType(path + '/.ipython/profile_profile_xfit/startup//ui/shell2.ui')):
    def __init__(self):
        super().__init__()
        self.setupUi(self)


# if __name__ == '__main__':
#     app = QtWidgets.QApplication([])
#     window = Shell()
#     window.show()
#     sys.exit(app.exec_())