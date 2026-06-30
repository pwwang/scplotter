# Shared argument definitions for spatial plot functions

This documentation block defines the full set of parameters shared
across all spatial plotting functions in scplotter:
[`SpatFeaturePlot`](https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.md),
[`SpatDimPlot`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md),
and the internal
[`SpatPlot()`](https://pwwang.github.io/scplotter/reference/SpatPlot.md)
workhorse. Each function inherits these parameter definitions via
`@inheritParams`. Individual function documentation pages may override
or supplement parameter descriptions when behavior differs between
functions.

## Arguments

- object:

  A Seurat object (with spatial data loaded via SeuratObject) or a
  Giotto object (created with Giotto). The spatial technology is
  auto-detected from the object's image class.

- fov:

  The name of the field of view (FOV) to plot. For Seurat FOV-based
  objects (Xenium, CosMx, etc.), defaults to
  [`SeuratObject::DefaultFOV()`](https://satijalab.github.io/seurat-object/reference/DefaultFOV.html).
  Not applicable to Visium or Slide-seq objects.

- boundaries:

  The name of the segmentation boundaries within the FOV to use for cell
  outlines. For Seurat FOV-based objects, defaults to
  [`SeuratObject::DefaultBoundary()`](https://satijalab.github.io/seurat-object/reference/Boundaries.html).
  Not applicable to Visium or Slide-seq objects.

- image:

  Controls the image/background layer. Possible values:

  - `NULL` — Default behavior: for Visium, the first image is used; for
    Giotto and FOV objects, no image is plotted.

  - A character string naming an image in the object — that specific
    image is plotted.

  - A color name (e.g., `"white"`, `"lightgray"`) — fills the background
    with a solid color rectangle.

  - `TRUE` — For Visium: uses the first image. For Giotto FOV: plots all
    non-overlapping images. For Seurat FOV: raises an error (no single
    default image).

  - `FALSE` — Disables the image layer entirely.

- masks:

  Logical. Whether to plot cell segmentation masks. Currently not
  supported — setting this to `TRUE` will produce an error for all
  object types.

- shapes:

  Controls the cell boundary (shapes) layer. Possible values:

  - `TRUE` — For Seurat FOV objects, uses `boundaries` as the shape
    boundaries. For other object types, uses the default boundaries. Not
    supported for Visium or Slide-seq objects.

  - A character string — The name of a specific set of boundaries within
    the FOV (Seurat) or spatial info name (Giotto) to use as shapes.

  - `FALSE` — Disables the shapes layer.

  Defaults to `TRUE` when `shapes_fill_by` is provided, and `NULL`
  otherwise.

- points:

  Logical. Whether to plot the points layer (cells, spots, or molecules
  as points on the spatial coordinates). Default: `TRUE`.

- ext:

  The spatial extent (bounding box) of the plot. If `NULL`, the extent
  is calculated automatically from the data or, when `crop = TRUE`, from
  the tissue coordinates. Can be a numeric vector in the format
  `c(xmin, xmax, ymin, ymax)` or a
  [`terra::SpatExtent`](https://rspatial.github.io/terra/reference/ext.html)
  object.

- crop:

  Logical. Whether to crop the plot to the extent of the tissue/spots.
  When `TRUE` (default), the plot is automatically zoomed to the data
  extent with optional `padding`. Analogous to the `crop` argument in
  [`Seurat::SpatialDimPlot()`](https://satijalab.org/seurat/reference/SpatialPlot.html).

- group_by:

  A metadata column name used to color the points. Must be a character
  or factor column in the object's metadata. For
  [`SpatDimPlot()`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md),
  if `NULL` and the object has FOV data with features, defaults to
  `"molecules"`; otherwise defaults to `"Identity"` (the active cluster
  identities). For
  [`SpatFeaturePlot()`](https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.md),
  `group_by` is ignored — use
  [`SpatDimPlot`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md)
  for categorical grouping. The special value `"molecules"` enables
  molecule-level visualization in FOV-based data.

- features:

  A character vector of feature names to visualize. For
  [`SpatFeaturePlot()`](https://pwwang.github.io/scplotter/reference/SpatFeaturePlot.md),
  each feature is plotted as a separate facet (or combined theme when a
  single feature is given), with expression values coloring the points.
  For
  [`SpatDimPlot()`](https://pwwang.github.io/scplotter/reference/SpatDimPlot.md),
  features are treated as molecule names to plot at single-molecule
  resolution (FOV data only). Can include gene names, metadata column
  names, or dimension reduction components.

- layer:

  The assay layer from which to extract feature expression values. For
  Seurat objects, one of `"data"` (default), `"scale.data"`, or
  `"counts"`. For Giotto objects, one of `"normalized"` (default),
  `"scaled"`, `"raw"`, `"counts"`, or `"custom"`.

- scale_factor:

  Internal use. The image scale factor extracted from the object, used
  to map between pixel and tissue coordinate spaces. Automatically
  determined from the object's image data.

- layers:

  A character vector specifying which layers to include and in what
  order. Possible values are `"image"`, `"masks"`, `"shapes"`, and
  `"points"`. Order matters — the first element is drawn first (bottom).
  Omit a layer name to disable it. Default:
  `c("image", "masks", "shapes", "points")` intersected with which
  layers are non-`FALSE`/non-`NULL`.

- flip_y:

  Logical. Whether to flip the y-axis. This is primarily for internal
  coordinate system alignment — Visium/Slide-seq objects default to
  `TRUE`, FOV-based objects to `FALSE`. In most cases you do not need to
  set this manually.

- padding:

  Numeric. Extra space added around the data extent when `crop = TRUE`.
  Expressed as a fraction of the data range (e.g., `0.05` adds 5\\
  defaults to `0`. For other object types, defaults to `0` when an image
  layer is present and `0.05` otherwise.

- image_scale:

  The image scale factor name to use, typically `"lowres"` or `"hires"`.
  Controls which resolution of the stored image is rendered. Analogous
  to the `image.scale` argument in
  [`Seurat::SpatialDimPlot()`](https://satijalab.org/seurat/reference/SpatialPlot.html).

- x:

  Internal use. The name of the x-coordinate column in the spatial data.
  Auto-detected based on the spatial technology (`"imagerow"` for Visium
  V1, varies for other types).

- y:

  Internal use. The name of the y-coordinate column in the spatial data.
  Auto-detected based on the spatial technology (`"imagecol"` for Visium
  V1, varies for other types).

- nmols:

  Integer. Maximum number of molecules to plot per feature when
  `group_by = "molecules"`. Analogous to the `nmols` argument in
  [`Seurat::ImageDimPlot()`](https://satijalab.org/seurat/reference/ImageDimPlot.html).
  Default: `1000`.

- shapes_fill_by:

  A column name in the metadata (or a feature/gene name) used to fill
  the cell boundary shapes. If a single color string is provided, all
  shapes are filled with that color. When set, `shapes` defaults to
  `TRUE`.

- graph:

  The name of a spatial network graph to overlay on the plot. Currently
  supported only for Giotto objects. Possible values:

  - `TRUE` — Use the default spatial network.

  - A character string — The graph name. If the name contains `":"`, the
    part before the colon is used as `spat_unit` and the part after as
    the graph name.

  - `NULL` (default) — No graph overlay.

  The graph data is retrieved via
  [`GiottoClass::getSpatialNetwork()`](https://giotto-suite.github.io/GiottoClass/reference/getSpatialNetwork.html).

- shape:

  Numeric. The point shape (ggplot2 shape aesthetic). Default: `16`
  (filled circle). See
  <https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html>
  for the full shape palette.

- legend.position:

  Character. Legend position. One of `"right"` (default), `"left"`,
  `"top"`, `"bottom"`, or `"none"`.

- legend.direction:

  Character. Legend direction. One of `"vertical"` (default) or
  `"horizontal"`.

- theme:

  A theme function or a character string naming one. Default:
  `"theme_box"`. Built-in aliases (usable without namespace):
  `"theme_box"`
  ([`plotthis::theme_box()`](https://pwwang.github.io/plotthis/reference/theme_box.html)),
  `"theme_this"`
  ([`plotthis::theme_this()`](https://pwwang.github.io/plotthis/reference/theme_this.html)),
  `"theme_blank"`
  ([`ggplot2::theme_void()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).
  Any ggplot2 theme can be used with its fully qualified name (e.g.,
  `"ggplot2::theme_bw"`).

- theme_args:

  A named list of additional arguments passed to the theme function.
  Default: [`list()`](https://rdrr.io/r/base/list.html).

- title:

  Character. Plot title. Default: `NULL` (no title).

- subtitle:

  Character. Plot subtitle. Default: `NULL`.

- xlab:

  Character. x-axis label. Default: `NULL`.

- ylab:

  Character. y-axis label. Default: `NULL`.

- facet_scales:

  Character. Whether facet scales are `"fixed"` (default), `"free"`,
  `"free_x"`, or `"free_y"`. Passed to
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  when multiple features are plotted.

- facet_nrow:

  Integer. Number of facet rows. Default: `NULL` (auto-calculated).

- facet_ncol:

  Integer. Number of facet columns. Default: `NULL` (auto-calculated).

- facet_byrow:

  Logical. Whether to fill facets by row. Default: `TRUE`.

- feat_type:

  Character. The feature type (modality) to query for expression values
  in Giotto objects. Common values: `"rna"` (default), `"dna"`,
  `"protein"`. Ignored for Seurat objects.

- use_overlap:

  Logical. For Giotto FOV objects, whether to use pre-computed
  polygon-feature overlap results (from
  [`GiottoClass::combineFeatureOverlapData()`](https://giotto-suite.github.io/GiottoClass/reference/combineFeatureOverlapData.html))
  instead of cell-level expression. Default: `FALSE`.

- shapes_feat_type:

  Character. The feature type to use when extracting metadata for shape
  filling in Giotto objects. Default: `"cell"`.

- shapes_alpha:

  Numeric. Transparency (alpha) value for the shapes layer, between 0
  and 1. When points are also plotted, defaults to `0.5` so points
  remain visible on top of shapes; otherwise defaults to `1`.

- spat_unit:

  Character. The spatial unit to query in a Giotto object (e.g.,
  `"cell"`, `"subcellular"`). Auto-detected if `NULL`. Ignored for
  Seurat objects.

- spat_loc_name:

  Character. The spatial locations name to query in a Giotto object.
  Auto-detected from available spatial locations if `NULL`. Ignored for
  Seurat objects.

- spat_enr_names:

  Character. Spatial enrichment results names in a Giotto object (for
  enrichment-based feature extraction). Ignored for Seurat objects.

- ...:

  Additional arguments passed to the underlying layer functions.
  Arguments are dispatched by prefix:

  `image_*`

  :   Arguments passed to
      [`plotthis::SpatImagePlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html)
      (e.g., `image_alpha`, `image_interpolation`).

  `masks_*`

  :   Arguments passed to
      [`plotthis::SpatMasksPlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html).

  `shapes_*`

  :   Arguments passed to
      [`plotthis::SpatShapesPlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html)
      (e.g., `shapes_color`, `shapes_linewidth`).

  `points_*`

  :   Arguments passed to
      [`plotthis::SpatPointsPlot()`](https://pwwang.github.io/plotthis/reference/spatialplots.html)
      (e.g., `points_size`, `points_alpha`).

  No prefix

  :   Arguments without a recognized prefix are treated as points
      arguments, but with lower priority than `points_*` arguments.

## Spatial data technologies

scplotter supports several spatial transcriptomics technologies through
both Seurat and Giotto objects:

- **Visium (10x Genomics)** — Spot-based spatial transcriptomics.
  Coordinates are stored with `imagerow`/`imagecol` axes. The tissue
  image is always available and plotted as the background layer.
  Supported via `SpatPlot.Seurat.Visium`.

- **Slide-seq** — Bead-based spatial transcriptomics with randomized
  barcode positions. Similar to Visium in coordinate handling but
  without a registered tissue image. Supported via
  `SpatPlot.Seurat.SlideSeq`.

- **FOV-based (Xenium, CosMx, MERSCOPE, Nanostring)** — Single-molecule
  resolution technologies with cell segmentation boundaries. Supports
  molecule-level visualization, cell boundary shapes, and registered
  images. Supported via `SpatPlot.Seurat.FOV` and `SpatPlot.giotto`.

## Layer system

Spatial plots in scplotter are composed of ordered layers. The `layers`
argument controls which layers are drawn and in what order:

1.  **"image"** — The tissue image or a colored background rectangle.
    Must be the first layer when included. Use `image = "colorname"` to
    fill the background with a solid color instead of an actual image.

2.  **"masks"** — Cell segmentation masks (not yet supported for Seurat
    or Giotto objects).

3.  **"shapes"** — Cell boundary polygons, filled by metadata or feature
    expression when `shapes_fill_by` is set.

4.  **"points"** — Individual points representing cells, spots, or
    molecules, colored by `group_by`, `features`, or `shapes_fill_by`.

Layer-specific styling is controlled via prefixed `...` arguments:
`image_*` for the image layer, `masks_*` for masks, `shapes_*` for
shapes, and `points_*` for points.

## Coordinate systems

Different spatial technologies use different coordinate conventions:

- **Visium/Slide-seq** — Coordinates are swapped (x and y are reversed)
  relative to the tissue image, and the y-axis is flipped
  (`flip_y = TRUE`) so the image origin aligns with the plot origin.

- **FOV-based (Seurat)** — Coordinates are in image space with no
  swapping needed. `flip_y = FALSE` by default.

- **FOV-based (Giotto)** — Coordinates are stored as `x` and `y`
  directly. `flip_y = FALSE` by default.

These are handled automatically; the `x`, `y`, and `flip_y` arguments
are primarily for internal use.
