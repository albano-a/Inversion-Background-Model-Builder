from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QTableWidgetItem
from PyQt5.QtCore import QDir
from PyQt5 import uic

MB = "src/gui/qt5/AdmaModelBuilder.ui"


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        uic.loadUi(MB, self)

        self.mbwAddGeometryButton.clicked.connect(self.add_geometry)
        self.mbwAddWellButton.clicked.connect(self.add_wells)
        self.mbwAddHorizonButton.clicked.connect(self.add_horizon)

        self.show()

    def add_geometry(self):
        pass

    def add_wells(self):
        wells = QFileDialog.getOpenFileNames(
            self, caption="Choose files", directory=QDir.currentPath()
        )[0]

        for well in wells:
            row_count = self.mbwWellTableWidget.rowCount()
            self.mbwWellTableWidget.insertRow(row_count)
            self.mbwWellTableWidget.setItem(row_count, 0, QTableWidgetItem(well))

    def add_horizon(self):
        horizons = QFileDialog.getOpenFileNames(
            self, caption="Choose files", directory=QDir.currentPath()
        )[0]

        for horizon in horizons:
            row_count = self.mbwHorizonTableWidget.rowCount()
            self.mbwHorizonTableWidget.insertRow(row_count)
            self.mbwHorizonTableWidget.setItem(row_count, 0, QTableWidgetItem(horizon))


if __name__ == "__main__":
    app = QApplication([])
    window = MainWindow()

    app.exec_()
