#!/usr/bin/env python3
"""Render the editable workflow SVG with the paper-wide DejaVu Sans typeface."""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


GROUP_STYLE = {
    "teal": ("#E1F5EE", "#0F6E56", {"th": "#085041", "ts": "#0F6E56"}),
    "coral": ("#FAECE7", "#993C1D", {"th": "#712B13", "ts": "#993C1D"}),
    "purple": ("#EEEDFE", "#534AB7", {"th": "#3C3489", "ts": "#534AB7"}),
    "gray": ("#F1EFE8", "#5F5E5A", {"th": "#2C2C2A", "ts": "#5F5E5A"}),
    "method": ("#E7ECEF", "#5A6B73", {"th": "#2C3E45", "ts": "#5A6B73", "tss": "#5A6B73"}),
    "sel": ("#FBF0CE", "#B8860B", {"th": "#6B4E00", "ts": "#8A6D0B"}),
}
TEXT_STYLE = {
    "t": (14.0, "#2C2C2A", "normal"),
    "ts": (12.0, "#5F5E5A", "normal"),
    "tss": (10.5, "#5F5E5A", "normal"),
    "th": (14.0, "#2C2C2A", "medium"),
}


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def number(value: str | None, default: float = 0.0) -> float:
    if value is None:
        return default
    return float(value.removesuffix("px"))


def inline_style(value: str | None) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in (value or "").split(";"):
        if ":" in item:
            key, val = item.split(":", 1)
            result[key.strip()] = val.strip()
    return result


def render(source: Path, output: Path) -> None:
    root = ET.parse(source).getroot()
    width = number(root.get("width"))
    height = number(root.get("height"))

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.sans-serif": ["DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(width / 96.0, height / 96.0), dpi=96, facecolor="white")
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.axis("off")

    def walk(element: ET.Element, group: str | None = None) -> None:
        tag = local_name(element.tag)
        current_group = element.get("class") if tag == "g" else group

        if tag == "rect":
            x, y = number(element.get("x")), number(element.get("y"))
            w, h = number(element.get("width")), number(element.get("height"))
            radius = number(element.get("rx"))
            if current_group in GROUP_STYLE:
                fill, edge, _ = GROUP_STYLE[current_group]
                lw = 1.4 if current_group == "sel" else 0.8
            else:
                fill = element.get("fill", "none")
                edge = element.get("stroke", "none")
                lw = number(element.get("stroke-width"), 0.8)
            patch = FancyBboxPatch(
                (x, height - y - h),
                w,
                h,
                boxstyle=f"round,pad=0,rounding_size={radius}",
                facecolor=fill,
                edgecolor=edge,
                linewidth=lw * 0.75,
            )
            ax.add_patch(patch)

        elif tag == "line":
            x1, y1 = number(element.get("x1")), height - number(element.get("y1"))
            x2, y2 = number(element.get("x2")), height - number(element.get("y2"))
            lw = number(element.get("stroke-width"), 1.0) * 0.75
            if element.get("marker-end"):
                ax.add_patch(
                    FancyArrowPatch(
                        (x1, y1),
                        (x2, y2),
                        arrowstyle="->",
                        mutation_scale=7,
                        linewidth=lw,
                        color="#888780",
                        shrinkA=0,
                        shrinkB=0,
                    )
                )
            else:
                ax.plot((x1, x2), (y1, y2), color="#888780", linewidth=lw)

        elif tag == "text":
            cls = element.get("class", "t")
            size_px, color, weight = TEXT_STYLE.get(cls, TEXT_STYLE["t"])
            overrides = inline_style(element.get("style"))
            if "font-size" in overrides:
                size_px = number(overrides["font-size"])
            if current_group in GROUP_STYLE:
                color = GROUP_STYLE[current_group][2].get(cls, color)
            anchor = element.get("text-anchor", "start")
            ax.text(
                number(element.get("x")),
                height - number(element.get("y")),
                "".join(element.itertext()),
                ha="center" if anchor == "middle" else "left",
                va="baseline",
                fontsize=size_px * 0.75,
                fontfamily="DejaVu Sans",
                fontweight=weight,
                color=color,
            )

        for child in element:
            if local_name(child.tag) not in {"style", "title", "desc", "defs", "marker", "path"}:
                walk(child, current_group)

    walk(root)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output,
        format="pdf",
        dpi=96,
        facecolor="white",
        edgecolor="none",
        metadata={"Creator": "ORT deterministic workflow renderer", "CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    render(args.source, args.output)


if __name__ == "__main__":
    main()
