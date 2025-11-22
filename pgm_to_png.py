import pathlib
import imageio.v2 as iio

def pgm_to_png_dir(folder):
    folder = pathlib.Path(folder)
    for pgm_path in folder.glob("*.pgm"):
        png_path = pgm_path.with_suffix(".png")
        img = iio.imread(pgm_path)          # lee PGM en array
        iio.imwrite(png_path, img)          # escribe PNG
        pgm_path.unlink()                   # elimina el PGM
        print(f"{pgm_path.name} -> {png_path.name}")

if __name__ == "__main__":
    pgm_to_png_dir("results")  # usa la carpeta actual; cámbialo si quieres otra
