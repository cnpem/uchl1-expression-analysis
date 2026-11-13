"""Figure 2 canvas widget: a Vega-Lite chart that reports clicked variants and a crop brush to Python.

`mo.ui.altair_chart` (marimo 0.25) does not report selections from layered charts, and the
lollipop canvas is layered (stems, heads, picked-variant rings). This anywidget renders the same
Altair spec with vega-embed and syncs two values back:

* `selected`: variant IDs; a click on a variant mark toggles its ID (drags that pan are ignored).
* `crop`: `[start, end]` in bp, set when a drag on the `crop_brush` interval ends.

The chart menu ("...") saves the figure as SVG or as PNG at 4x resolution, for papers and posters.

Nothing here is specific to one gene.
"""

from __future__ import annotations

import anywidget
import traitlets

from viz.workbench_helpers import CROP_BRUSH_NAME

VEGA_EMBED_URL = "https://esm.sh/vega-embed@7?deps=vega@6,vega-lite@6.4.1"

_ESM = """
import vegaEmbed from "%(url)s";

function render({ model, el }) {
  let view = null;
  let down = null;
  let pendingCrop = null;
  let timer = null;

  async function draw() {
    if (view) { view.finalize(); view = null; }
    el.innerHTML = "";
    const host = document.createElement("div");
    host.style.width = "100%%";
    el.appendChild(host);
    const result = await vegaEmbed(host, model.get("spec"), {
      renderer: "canvas",
      // "..." menu: Save as SVG (vector, for papers) or PNG at 4x resolution (posters)
      actions: { export: true, source: false, compiled: false, editor: false },
      scaleFactor: { png: 4, svg: 1 },
      downloadFileName: model.get("file_name"),
    });
    view = result.view;
    view.addEventListener("click", (event, item) => {
      if (down && Math.hypot(event.clientX - down[0], event.clientY - down[1]) > 4) return;
      const id = item && item.datum ? item.datum.variant_id : undefined;
      if (id === undefined || id === null) return;
      const current = [...model.get("selected")];
      const i = current.indexOf(id);
      if (i >= 0) current.splice(i, 1); else current.push(id);
      model.set("selected", current);
      model.save_changes();
    });
    try {
      view.addSignalListener("%(brush)s", (name, value) => {
        const bounds = value && Object.values(value)[0];
        pendingCrop = Array.isArray(bounds) && bounds.length === 2 ? bounds : null;
        clearTimeout(timer);
        if (pendingCrop) timer = setTimeout(commitCrop, 700);
      });
    } catch (e) { /* no crop brush in this spec (pan mode) */ }
  }

  function commitCrop() {
    clearTimeout(timer);
    if (!pendingCrop) return;
    const [a, b] = pendingCrop;
    pendingCrop = null;
    const lo = Math.round(Math.min(a, b)), hi = Math.round(Math.max(a, b));
    if (hi - lo < model.get("min_crop_bp")) return;
    model.set("crop", [lo, hi]);
    model.save_changes();
  }

  // A crop is committed when the drag ends, or 0.7 s after the brush stops changing
  el.addEventListener("pointerdown", (e) => { down = [e.clientX, e.clientY]; pendingCrop = null; clearTimeout(timer); });
  const onUp = () => { if (pendingCrop) commitCrop(); };
  window.addEventListener("pointerup", onUp);
  model.on("change:spec", draw);
  draw();
  return () => { window.removeEventListener("pointerup", onUp); clearTimeout(timer); if (view) view.finalize(); };
}

export default { render };
""" % {"url": VEGA_EMBED_URL, "brush": CROP_BRUSH_NAME}


class VariantCanvas(anywidget.AnyWidget):
    _esm = _ESM
    spec = traitlets.Dict().tag(sync=True)
    selected = traitlets.List(traitlets.Unicode()).tag(sync=True)
    crop = traitlets.List(traitlets.Int()).tag(sync=True)
    min_crop_bp = traitlets.Int(5).tag(sync=True)
    file_name = traitlets.Unicode("figure2").tag(sync=True)


def chart_spec(chart) -> dict:
    """Vega-Lite dict with the data inlined (independent of marimo's data transformer)."""
    import altair as alt

    with alt.data_transformers.enable("default", max_rows=None):
        return chart.to_dict()


__all__ = ["CROP_BRUSH_NAME", "VariantCanvas", "chart_spec"]
