import io
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont


SHAPE = (1024, 1024, 1024)
OFFSET = 478
VOXEL_MM = 0.0719999969
LOW, HIGH = 0.0, 9000.0
CROP = 512
PANEL = 780
HALF_SLAB = 2


def get_font(size, bold=False):
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    for root in ("/usr/share/fonts/truetype/dejavu", "/usr/share/fonts/dejavu"):
        path = os.path.join(root, name)
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def crop(a, cx, cy):
    h = CROP // 2
    x0, x1 = int(round(cx)) - h, int(round(cx)) + h
    y0, y1 = int(round(cy)) - h, int(round(cy)) + h
    out = np.zeros((CROP, CROP), dtype=np.float32)
    sx0, sx1 = max(0, x0), min(a.shape[1], x1)
    sy0, sy1 = max(0, y0), min(a.shape[0], y1)
    out[sy0-y0:sy1-y0, sx0-x0:sx1-x0] = a[sy0:sy1, sx0:sx1]
    return out


def get_sections(path, x, y, z):
    v = np.memmap(path, dtype="<u2", mode="r", offset=OFFSET, shape=SHAPE)
    iz, iy = int(round(z)), int(round(y))
    plane1 = np.asarray(v[iz-HALF_SLAB:iz+HALF_SLAB+1, :, :], dtype=np.float32).mean(axis=0)
    plane2 = np.asarray(v[:, iy-HALF_SLAB:iy+HALF_SLAB+1, :], dtype=np.float32).mean(axis=1)
    return crop(plane1, x, y), crop(plane2, x, z)


def render(a, letter):
    u8 = np.clip((a - LOW) * 255.0 / (HIGH - LOW), 0, 255).astype(np.uint8)
    im = Image.fromarray(u8).resize((PANEL, PANEL), Image.Resampling.LANCZOS).convert("RGB")
    d = ImageDraw.Draw(im)
    d.rectangle((18, 18, 73, 73), fill=(255, 255, 255))
    d.text((30, 20), letter, fill=(15, 15, 15), font=get_font(38, True))
    bar_px = round((5.0 / (CROP * VOXEL_MM)) * PANEL)
    x0, y0 = 34, PANEL - 38
    d.line((x0, y0, x0 + bar_px, y0), fill="white", width=8)
    d.text((x0, y0 - 41), "5 mm", fill="white", font=get_font(26, True))
    return im


def main():
    specs = [
        ("ORT15 (left petrous)", sys.argv[1], *map(float, sys.argv[2:5])),
        ("ORT16 (right petrous)", sys.argv[5], *map(float, sys.argv[6:9])),
    ]
    data = [(label, get_sections(path, x, y, z)) for label, path, x, y, z in specs]
    left, top, gx, gy = 40, 100, 46, 68
    width = left + 2 * PANEL + gx + 40
    height = top + 2 * PANEL + gy + 55
    out = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(out)
    letters = (("A", "C"), ("B", "D"))
    for col, (label, sections) in enumerate(data):
        x0 = left + col * (PANEL + gx)
        hf = get_font(34, True)
        bb = d.textbbox((0, 0), label, font=hf)
        d.text((x0 + (PANEL - (bb[2] - bb[0])) / 2, 38), label, fill=(20, 20, 20), font=hf)
        for row, section in enumerate(sections):
            y0 = top + row * (PANEL + gy)
            out.paste(render(section, letters[col][row]), (x0, y0))
    buf = io.BytesIO()
    out.save(buf, format="PNG", optimize=True, dpi=(300, 300))
    sys.stdout.buffer.write(buf.getvalue())


if __name__ == "__main__":
    main()
