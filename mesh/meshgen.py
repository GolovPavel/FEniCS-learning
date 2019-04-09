import sys
from typing import List, Dict

import h5py
import numpy as np
import pygmsh
import meshio
from PIL import Image
from tqdm import tqdm

from config import Config
from logger import SystemLogger


class Cylinder:
    u"""Вспомогательный класс, содержащий основные параметры цилиндра"""

    def __init__(self, radius, height, mesh_size):
        self.radius = radius
        self.height = height
        self.mesh_size = mesh_size


class GmshCylinder:
    u"""Вспомогательный класс, содержащий основные элементы цилиндра в gmsh"""

    def __init__(self, points, circles, circle_ll, lines, surface_ll, surfaces, plane_surfaces, surface_loop, volume):
        self.points = points
        self.circles = circles
        self.circles_ll = circle_ll
        self.lines = lines
        self.surface_ll = surface_ll
        self.surfaces = surfaces
        self.plane_surfaces = plane_surfaces
        self.surface_loop = surface_loop
        self.volume = volume


class MeshGenerator:
    u"""Класс, отвечающий за создание расчетной сетки"""

    u"""Описание переменных экземляров класса"""
    gmsh_cylinders: np.array  # Массив с gmsh элементами цилиндров
    geometry: pygmsh.built_in.Geometry  # Объект, определяющий gmsh геометрию
    points_2d: np.array  # Массив узлов 2d сетки
    points_3d: np.array  # Массив узлов 3d сетки
    cells_2d: Dict  # Ячейки 2d сетки
    cells_3d: Dict  # Ячейки 3d сетки
    img_pixels: object  # Объект с пикселями изображения
    pix_width: int  # Ширина изображения в пикселях
    pix_height: int  # Высота изображения в пикселях
    img_center: List  # Координаты центра изображения (индесы img_pixels)
    radius_pix: int  # Радиус изображения в пикселях
    radius: int  # Радиус рассматриваемой зоны в метрах
    scale: float  # Расстояние, приходящееся на половину пикселя в метрах
    phys_coords: np.array  # Массив физ. координат метности
    phys_types: np.array  # Массив типов местности для каждой физической координаты (phys_coords)
    k: np.array # Массив степенных коэффициентов для определения скорости ветра в зависимости от высоты и типа подстилающей поверхности
    marked_points_2d: np.array  # Массив промаркированных координат нижней поверхности сетки
    marked_points_3d: np.array  # Массив промаркированных координат 3d сетки

    def __init__(self, config_file):
        # Инициализация конфигурации
        self.cfg = Config(config_file)
        self.cfg.init_mesh()
        self.cfg.init_pic()
        self.cfg.init_soft()
        self.cfg.init_koef()

        # Инициализация логгирования
        self.logger = SystemLogger(config_file)

        # Инициализация внутренних компонентов
        self.__init_cylinders()
        self.__init_image()

    def __init_cylinders(self):
        self.cylinders = []
        for i in range(len(self.cfg.radius_list)):
            cylinder = Cylinder(
                self.cfg.radius_list[i],
                self.cfg.h,
                self.cfg.step[i]
            )
            self.cylinders.append(cylinder)

    def __init_image(self):
        u"""Инициализация параметров изображения"""
        rgb_img = Image.open(self.cfg.image_path).convert("RGB")
        self.img_pixels = rgb_img.load()
        self.pix_width = rgb_img.width
        self.pix_height = rgb_img.height
        self.img_center = [int(self.pix_width / 2), int(self.pix_height / 2)]
        self.radius_pix = int(min([self.pix_width, self.pix_height]) / 2)
        self.radius = self.cfg.radius_list[-1]
        self.scale = self.radius / (2 * self.radius_pix)

    def generate_and_save(self):
        self.logger.print_log("Начинаем процесс генерации расчетной сетки...", "INFO")
        self.__generate_mesh()
        self.__save_xml_mesh()
        self.__convert_image()
        self.__mark_mesh()
        self.__save_marked_mesh()
        self.logger.print_log("Процесс генерации расчетной сетки успешно завершен!", "INFO")
        self.logger.print_log("Путь к файлу расчетной сетки: %s" % self.cfg.h5_mesh_path, "INFO")

    def __generate_mesh(self):
        u"""Создания сетки на основе построенной gmsh геометрии.
        Создает 2 массива: координаты узлов сетки в 3d и координаты узлов нижней поверхности в 2d геометрии."""
        self.logger.print_log("Создаем сетку...", "INFO")
        self.__create_geometry()
        self.points_3d, self.cells_3d, _, _, _ = pygmsh.generate_mesh(self.geometry, dim=3, gmsh_path=self.cfg.gmsh_path)
        self.__mark_bottom_surface()  # Помечаем нижнюю поверхность геометрии для создания сетки только на ней
        self.points_2d, self.cells_2d, _, _, _ = pygmsh.generate_mesh(self.geometry, dim=2, gmsh_path=self.cfg.gmsh_path)
        self.logger.print_log("Сетка создана успешно.", "INFO")

    def __save_xml_mesh(self):
        u"""Сохранения расчетной сетки в dolfin-xml формат"""
        self.logger.print_log("Сохраняем сетку в xml формат...", "INFO")
        meshio.write_points_cells(self.cfg.xml_2d_path, self.points_2d, self.cells_2d)
        meshio.write_points_cells(self.cfg.xml_3d_path, self.points_3d, self.cells_3d)
        self.logger.print_log("Сетка успешно сохранена в xml формат", "INFO")

    def __convert_image(self):
        u"""Преобразование пикселей изображения в физические координаты с соответствующим типом местности"""
        self.logger.print_log("Преобразуем координаты пикселей в физические координаты...", "INFO")
        self.phys_coords = []
        self.phys_types = []
        for y_pix in range(self.pix_height):
            # Переводим координату y в систему координат с началом отсчета в центре изображения
            y_coord = self.radius / self.radius_pix * (self.img_center[1] - y_pix)
            y_coord = y_coord - self.scale if y_coord > 0 else y_coord + self.scale
            for x_pix in range(self.pix_width):
                # Обрезаем изображение по радиусу в пикселях
                if (x_pix - self.img_center[0]) ** 2 + (y_pix - self.img_center[1]) ** 2 <= self.radius_pix ** 2:
                    # Переводим координату x в систему координат с началом отсчета в центре изображения
                    x_coord = self.radius / self.radius_pix * (x_pix - self.img_center[0])
                    x_coord = x_coord - self.scale if x_coord > 0 else x_coord + self.scale
                    pixel_color = self.img_pixels[x_pix, y_pix]
                    pixel_type = self.__determine_pixel_type(pixel_color)
                    self.phys_coords.append((x_coord, y_coord))
                    self.phys_types.append(pixel_type)
        self.phys_coords = np.array(self.phys_coords)
        self.phys_types = np.array(self.phys_types)
        self.logger.print_log("Координаты пикселей успешно преобразованы в физические координаты", "INFO")

    def __mark_mesh(self):
        u"""Маркировка узлов 2d и 3d сеток по типу местности и высоте относительно подстилающей поверхности"""
        self.logger.print_log("Маркируем сетку по типу местности...", "INFO")
        self.marked_points_2d = np.zeros(self.points_2d.shape[0],
                                         dtype=[('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('surf_type', 'S10'), ('z0_position', 'S10')])
        self.marked_points_3d = np.zeros(self.points_3d.shape[0],
                                         dtype=[('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('surf_type', 'S10'), ('z0_position', 'S10'), ('k', 'f8')])
        self.__mark_2d_mesh()
        self.__mark_3d_mesh()
        self.logger.print_log("Маркировка сетки по типу местности успешно завершена.", "INFO")

    def __save_marked_mesh(self):
        u"""Запись координат сетки с индикаторами местиности в h5 файл"""
        with h5py.File(self.cfg.h5_mesh_path, 'w') as f:
            cols = [
                ("X", np.float64),
                ("Y", np.float64),
                ("Z", np.float64),
                ("SURF_TYPE", "S10"),
                ("Z0_POSITION", np.int8),
                ("K", np.float64),
            ]
            mesh_2d = np.zeros(len(self.marked_points_2d), dtype=cols)
            mesh_3d = np.zeros(len(self.marked_points_3d), dtype=cols)

            mesh_2d["X"][:] = [coord[0] for coord in self.marked_points_2d]
            mesh_2d["Y"][:] = [coord[1] for coord in self.marked_points_2d]
            mesh_2d["Z"][:] = [coord[2] for coord in self.marked_points_2d]
            mesh_2d["SURF_TYPE"][:] = [coord[3] for coord in self.marked_points_2d]
            mesh_2d["Z0_POSITION"][:] = [coord[4] for coord in self.marked_points_2d]
            mesh_2d["K"][:] = 0.0 # Не будет использоваться

            mesh_3d["X"][:] = [coord[0] for coord in self.marked_points_3d]
            mesh_3d["Y"][:] = [coord[1] for coord in self.marked_points_3d]
            mesh_3d["Z"][:] = [coord[2] for coord in self.marked_points_3d]
            mesh_3d["SURF_TYPE"][:] = [coord[3] for coord in self.marked_points_3d]
            mesh_3d["Z0_POSITION"][:] = [coord[4] for coord in self.marked_points_3d]
            mesh_3d["K"][:] = [coord[5] for coord in self.marked_points_3d]

            f.create_dataset(name="2D_COORDS", data=mesh_2d)
            f.create_dataset(name="3D_COORDS", data=mesh_3d)

    def __create_geometry(self):
        u"""Построения gmsh геометрии, состоящей из циклидров (массив self.cylinders)"""
        self.geometry = pygmsh.built_in.Geometry()
        self.gmsh_cylinders = []
        for cylinder in self.cylinders:
            self.gmsh_cylinders.append(
                self.__add_cylinder(cylinder, self.gmsh_cylinders[-1] if len(self.gmsh_cylinders) != 0 else None)
            )

    def __mark_bottom_surface(self):
        u"""Маркировка нижней поверхности геометрии"""
        bottom_surfaces = [cylinder.plane_surfaces[0] for cylinder in self.gmsh_cylinders]
        center_point = self.gmsh_cylinders[0].points[0]
        self.geometry.add_physical_surface(bottom_surfaces)
        self.geometry.add_physical_point(center_point)

    def __determine_pixel_type(self, pixel_color):
        u"""Определение типа местности по цвету пикселя"""
        if pixel_color == self.cfg.rgb["forest"]:
            return self.cfg.marks["forest"]
        elif pixel_color == self.cfg.rgb["water"]:
            return self.cfg.marks["water"]
        elif pixel_color == self.cfg.rgb["town"]:
            return self.cfg.marks["town"]
        else:
            return self.cfg.marks["field"]

    def __mark_2d_mesh(self):
        u"""Маркировка узлов нижней поверхности 2d сетки"""
        self.logger.print_log("Маркируем 2d сетку по типу местности...", "INFO")
        for i, mesh_point in enumerate(tqdm(self.points_2d)):
            delta_x = mesh_point[0] - self.phys_coords.T[0]
            delta_y = mesh_point[1] - self.phys_coords.T[1]
            distance = np.sqrt(np.power(delta_x, 2) + np.power(delta_y, 2))
            surface_type = self.phys_types[distance.argmin()]
            z0_position = 1
            self.marked_points_2d[i] = (*mesh_point, surface_type, z0_position)

    def __mark_3d_mesh(self):
        u"""Маркировка узлов нижней поверхности 3d сетки"""
        self.logger.print_log("Маркируем 3d сетку по типу местности...", "INFO")
        for i, mesh_point in enumerate(tqdm(self.points_3d)):
            delta_x = mesh_point[0] - self.phys_coords.T[0]
            delta_y = mesh_point[1] - self.phys_coords.T[1]
            distance = np.sqrt(np.power(delta_x, 2) + np.power(delta_y, 2))
            surface_type = self.phys_types[distance.argmin()]
            z0_position = 0 if mesh_point[2] > self.cfg.z0[surface_type] else 1
            k = self.cfg.k[surface_type]
            self.marked_points_3d[i] = (*mesh_point, surface_type, z0_position, k)

    def __add_cylinder(self, cylinder, prev_cylinder):
        u"""Добавление в текущую gmsh геометрию нового цилиндра
        Если создаваемый цилиндр внутренний, то prev_cylinder должен быть None, иначе в prev_cylinder передается
        объект класса GmshElements, содержащий элементы внутреннего цилиндра"""
        base_coords = [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, cylinder.height]),
        ]

        coords = [
            np.array([cylinder.radius, 0.0, 0.0]),
            np.array([0.0, cylinder.radius, 0.0]),
            np.array([-cylinder.radius, 0.0, 0.0]),
            np.array([0.0, -cylinder.radius, 0.0]),
            np.array([cylinder.radius, 0.0, cylinder.height]),
            np.array([0.0, cylinder.radius, cylinder.height]),
            np.array([-cylinder.radius, 0.0, cylinder.height]),
            np.array([0.0, -cylinder.radius, cylinder.height])
        ]

        circle_points_idx = [
            np.array([2, 0, 3]),
            np.array([3, 0, 4]),
            np.array([4, 0, 5]),
            np.array([5, 0, 2]),

            np.array([6, 1, 7]),
            np.array([7, 1, 8]),
            np.array([8, 1, 9]),
            np.array([9, 1, 6])
        ]

        circles_line_loop_idx = [
            np.array([0, 1, 2, 3]),
            np.array([4, 5, 6, 7]),
        ]

        line_points_idx = [
            np.array([6, 2]),
            np.array([7, 3]),
            np.array([8, 4]),
            np.array([9, 5]),
        ]

        surfaces_line_loop_idx = [
            np.array([0, 0, 1, 4]),
            np.array([1, 1, 2, 5]),
            np.array([2, 2, 3, 6]),
            np.array([3, 3, 0, 7])
        ]

        points = self.__create_points(base_coords, coords, cylinder.mesh_size, prev_cylinder)
        circles = self.__create_circles(points, circle_points_idx)
        circle_ll = self.__create_circle_line_loops(circles, circles_line_loop_idx)
        lines = self.__create_lines(points, line_points_idx)
        surface_ll = self.__create_surface_line_loops(lines, circles, surfaces_line_loop_idx)
        surfaces = self.__create_surfaces(surface_ll)
        plane_surfaces = self.__create_plane_surfaces(circle_ll, prev_cylinder)
        surface_loop = self.__create_surface_loop(surfaces, plane_surfaces, prev_cylinder)
        volume = self.__create_volume(surface_loop)

        return GmshCylinder(points, circles, circle_ll, lines, surface_ll,
                            surfaces, plane_surfaces, surface_loop, volume)

    def __create_points(self, base_coords, point_coords, lcar, prev_cylinder):
        u"""Создание опорных точек в gmsh геометрии"""
        points = [self.geometry.add_point(point, lcar=lcar) for point in point_coords]
        if prev_cylinder is not None:
            points = [prev_cylinder.points[0], prev_cylinder.points[1]] + points
        else:
            points = [self.geometry.add_point(point, lcar=0.0) for point in base_coords] + points
        return points

    def __create_circles(self, points, circle_points_idx):
        u"""Создание дуг окружностей в gmsh геометрии"""
        circles = [
            self.geometry.add_circle_arc(points[idx[0]], points[idx[1]], points[idx[2]])
            for idx in circle_points_idx
        ]
        return circles

    def __create_circle_line_loops(self, circles, circle_ll_idx):
        u"""Создание окружностей в gmsh геометрии"""
        circle_ll = [
            self.geometry.add_line_loop([circles[idx[0]], circles[idx[1]], circles[idx[2]], circles[idx[3]]])
            for idx in circle_ll_idx
        ]
        return circle_ll

    def __create_lines(self, points, line_points_idx):
        u"""Создание линий в gmsh геометрии"""
        lines = [self.geometry.add_line(points[idx[0]], points[idx[1]]) for idx in line_points_idx]
        return lines

    def __create_surface_line_loops(self, lines, circles, surfaces_line_loop_idx):
        u"""Объединение линий боковых поверхностей цилиндра в gmsh геометрии"""
        surface_line_loops = [
            self.geometry.add_line_loop([lines[idx[0]], circles[idx[1]], -lines[idx[2]], -circles[idx[3]]])
            for idx in surfaces_line_loop_idx
        ]
        return surface_line_loops

    def __create_surfaces(self, surface_line_loops):
        u"""Создание боковых поверхностей цилиндра в gmsh геометрии"""
        surfaces = [self.geometry.add_surface(surface) for surface in surface_line_loops]
        return surfaces

    def __create_plane_surfaces(self, circle_line_loops, prev_cylinder):
        u"""Создание поверхностей оснований цилиндра в gmsh геометрии"""
        if prev_cylinder is None:
            plane_surfaces = [self.geometry.add_plane_surface(circle_loop) for circle_loop in circle_line_loops]
        else:
            plane_surfaces = [
                self.geometry.add_plane_surface(circle_loop, [prev_cylinder.circles_ll[i]])
                for i, circle_loop in enumerate(circle_line_loops)
            ]
        return plane_surfaces

    def __create_surface_loop(self, surfaces, plane_surfaces, prev_cylinder):
        u"""Объединение всех поверхностей цилиндра между собой в gmsh геометрии"""
        if prev_cylinder is None:
            surface_loop = self.geometry.add_surface_loop([plane_surfaces[0], *surfaces, plane_surfaces[1]])
        else:
            surface_loop = self.geometry.add_surface_loop(
                [plane_surfaces[0], *surfaces, plane_surfaces[1], *prev_cylinder.surfaces]
            )
        return surface_loop

    def __create_volume(self, surface_loop):
        u"""Создание объема цилиндра в gmsh геометрии"""
        volume = self.geometry.add_volume(surface_loop)
        return volume


if __name__ == '__main__':
    config_file_name = sys.argv[1]
    mesh_generator = MeshGenerator(config_file_name)
    mesh_generator.generate_and_save()
