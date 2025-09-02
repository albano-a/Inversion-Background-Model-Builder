from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QTableWidgetItem
from PyQt5.QtCore import QDir
from PyQt5 import uic
from core import bmtools

MB = "src/gui/qt5/AdmaModelBuilder.ui"
_WELL_POSITIONS_SEISMIC = {
    "1-RJS-342-RJS": {"X": 1118, "Y": 1410},
    "3-RJS-355-RJS": {"X": 1387, "Y": 1347},
    "3-RJS-360A-RJS": {"X": 1220, "Y": 1703},
    "3-RJS-510A-RJS": {"X": 940, "Y": 2376},
    "4-ABL-30A-RJS": {"X": 1318, "Y": 2523},
    "4-RJS-367-RJS": {"X": 906, "Y": 2183},
    "4-RJS-477A-RJS": {"X": 889, "Y": 1584},
    "6-ABL-1-RJS": {"X": 1304, "Y": 2088},
    "9-AB-67-RJS": {"X": 1107, "Y": 1493},
    "9-AB-71-RJS": {"X": 1054, "Y": 1500},
    "9-ABL-2-RJS": {"X": 1359, "Y": 1694},
    "9-ABL-3B-RJS": {"X": 1017, "Y": 1778},
    "9-ABL-5-RJS": {"X": 1347, "Y": 1816},
    "9-ABL-6A-RJS": {"X": 1466, "Y": 1823},
    "9-ABL-7-RJS": {"X": 1190, "Y": 1780},
    "9-ABL-9D-RJS": {"X": 1071, "Y": 2031},
}


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        uic.loadUi(MB, self)

        self.mbwAddGeometryButton.clicked.connect(self.add_geometry)
        self.mbwAddWellButton.clicked.connect(self.add_wells)
        self.mbwAddHorizonButton.clicked.connect(self.add_horizon)
        self.mbwConfirmButton.clicked.connect(self.prepare_horizons)

        self.show()

    def add_geometry(self):
        geometries = QFileDialog.getOpenFileNames(
            self, caption="Choose files", directory=QDir.currentPath()
        )[0]

        for geometry in geometries:
            row_count = self.mbwGeometryTableWidget.rowCount()
            self.mbwGeometryTableWidget.insertRow(row_count)
            self.mbwGeometryTableWidget.setItem(
                row_count, 0, QTableWidgetItem(geometry)
            )

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

    def prepare_wells(self):
        positions = []
        for row in range(self.mbwWellTableWidget.rowCount()):
            well_path_item = self.mbwWellTableWidget.item(row, 0)
            if well_path_item:
                well_path = well_path_item.text()
                # Extract well name from path (assuming it's the filename without extension)
                well_name = well_path.split("/")[-1].split(".")[0]

                if well_name in _WELL_POSITIONS_SEISMIC:
                    x = _WELL_POSITIONS_SEISMIC[well_name]["X"]
                    y = _WELL_POSITIONS_SEISMIC[well_name]["Y"]
                    positions.append([x, y])

        return positions

    def prepare_horizons(self):
        horizons = []
        for row in range(self.mbwHorizonTableWidget.rowCount()):
            horizon_path_item = self.mbwHorizonTableWidget.item(row, 0)
            if horizon_path_item:
                horizon_path = horizon_path_item.text()
                grid = bmtools.load_and_grid(horizon_path)
                horizons.append(grid)

        print(horizons)
        return horizons

    def run_background_model(self, horizons, samples, velocity_flat):
        well_p = self.prepare_wells()
        well_logs = bmtools.stacked_wells_matrix("data/Wells", sample=4.0)

        model = bmtools.run_background_modelling(
            horizons=horizons,
            log_well=well_logs,
            well_positions=well_p,
            samples=samples,
            velocity_flat=velocity_flat,
            smooth_model=True,
        )


if __name__ == "__main__":
    app = QApplication([])
    window = MainWindow()

    app.exec_()
