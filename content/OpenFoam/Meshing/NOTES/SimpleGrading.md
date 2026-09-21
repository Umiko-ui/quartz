# blockMesh grading trace: `simpleGrading (4 1 1)` (OpenFOAM v2412)

All code below is copied from the files pasted from `/usr/lib/openfoam/openfoam2412/src/mesh/blockMesh/`. Shortcut: `export BM=$FOAM_SRC/mesh/blockMesh`

Legend: **[SEEN]** = exact code from your tree. **[NOT SEEN]** = not yet pasted, so nothing below is verified.

## 0. Pipeline overview

```
simpleGrading (4 1 1)
 1  blockMesh ctor -> createTopology()                       blockMesh/blockMesh.C
 2  createTopology(): blocks via block::iNew                 blockMesh/blockMeshTopology.C
 3  blockDescriptor ctor: hex, (20 20 20), skip keyword,
    List<gradingDescriptors>(is)                             blockDescriptor/blockDescriptor.C
 4  gradingDescriptors operator>>: 4 -> [{1,1,4}]            gradingDescriptor/gradingDescriptors.C
 5  assignGradings(): expand_[0..3]=[4], [4..7]=[1], [8..11]=[1]
                                                             blockDescriptor/blockDescriptor.C
    -- points() -> createPoints()  (trigger chain: see section 6) --
 6  block::createPoints(): edgesPointsWeights(p, w)          blocks/block/blockCreate.C
 7  calcEdgePointsWeights(): lineDivide(edge, N, expand)     blockDescriptor/blockDescriptorEdges.C
 8  lineDivide ctor: lambda_j = (1-r^j)/(1-r^N)              blockEdges/lineDivide/lineDivide.C
 9  points_ = cedge.position(lambda)                         blockEdges/lineEdge/lineEdge.C [NOT SEEN]
10  createPoints(): w[e][i]=lambda -> blended -> points_     blocks/block/blockCreate.C
```

Key idea: the dictionary only stores a ratio. Coordinates are computed later, in `lineDivide`.

---

## 1. Entry: `blockMesh` constructor [SEEN]

**File:** `blockMesh/blockMesh.C`

```cpp
Foam::blockMesh::blockMesh
(
    const IOdictionary& dict,
    const word& regionName,
    mergeStrategy strategy,
    int verbosity
)
:
    meshDict_(dict),
    ...
    blockVertices_
    (
        meshDict_.lookup("vertices"),
        blockVertex::iNew(meshDict_, geometry_)
    ),
    vertices_(Foam::vertices(blockVertices_)),
    prescaling_(vector::uniform(1)),
    scaling_(vector::uniform(1)),
    transform_(),
    topologyPtr_(createTopology(meshDict_, regionName))   // <-- hand-off
{
    ...
    if (mergeStrategy_ == mergeStrategy::MERGE_POINTS)
    {
        calcGeometricalMerge();
    }
    else
    {
        calcTopologicalMerge();
    }
}
```

Lazy accessors in the same file:

```cpp
const Foam::pointField& Foam::blockMesh::points() const
{
    if (points_.empty())
    {
        createPoints();
    }

    return points_;
}
```

`scale` (formerly `convertToMeters`), `prescale` and `transform` are applied after the points exist (`inplacePointTransforms`), so they never affect grading.

---

## 2. Blocks are created from the dict [SEEN]

**File:** `blockMesh/blockMeshTopology.C`, function `blockMesh::createTopology`

```cpp
    // Create the blocks
    if (verbose_)
    {
        Info<< "Creating topology blocks" << endl;
    }
    {
        blockList blocks
        (
            meshDescription.lookup("blocks"),
            block::iNew(meshDescription, vertices_, edges_, faces_)
        );

        transfer(blocks);
    }
```

`block::iNew` is in `blocks/block/block.H` [NOT SEEN]. It forwards the stream to the `block` constructor, which builds a `blockDescriptor`. Confirm with:

```bash
grep -n "iNew" $BM/blocks/block/block.H $BM/blocks/block/block.C
```

Edges are read just before this in the same function:

```cpp
        blockEdgeList edges
        (
            meshDescription.lookup("edges"),
            blockEdge::iNew(meshDescription, geometry_, vertices_)
        );
```

---

## 3. Parsing one block [SEEN]

**File:** `blockDescriptor/blockDescriptor.C`, constructor taking `Istream& is`

For `hex (0 1 2 3 4 5 6 7) (20 20 20) simpleGrading (4 1 1)`:

```cpp
    // Read cell model and list of vertices (potentially with variables)
    word model(is);                                   // "hex"
    blockShape_ = cellShape
    (
        model,
        blockMeshTools::read<label>
        (
            is,
            dict.subOrEmptyDict("namedVertices")
        )
    );                                                // (0 1 2 3 4 5 6 7)

    // Examine next token
    token t(is);

    // Optional zone name
    if (t.isWord())
    {
        zoneName_ = t.wordToken();

        // Examine next token
        is >> t;
    }
    is.putBack(t);

    if (t.isPunctuation())
    {
        // New-style: read a list of 3 values
        if (t.pToken() == token::BEGIN_LIST)
        {
            is >> ijkMesh::sizes();                   // (20 20 20)
        }
        ...
    }
    ...

    is >> t;
    if (!t.isWord())
    {
        is.putBack(t);
    }                                                 // "simpleGrading" is read and DISCARDED

    List<gradingDescriptors> expand(is);              // reads (4 1 1)

    if (!assignGradings(expand))
    {
        FatalErrorInFunction
            << "Unknown definition of expansion ratios: " << expand
            << exit(FatalError);
    }

    check(is);

    findCurvedFaces(blockIndex);
```

Consequences:

- `simpleGrading` and `edgeGrading` are interchangeable keywords. The **number of values** decides the meaning.
- N = 20 is stored in `ijkMesh::sizes()`. It is used later as `sizes()[edgei/4]`.

---

## 4. One grading number becomes a section [SEEN]

### 4a. One section: `gradingDescriptor`

**File:** `gradingDescriptor/gradingDescriptor.H` (data)

```cpp
        scalar blockFraction_;
        scalar nDivFraction_;
        scalar expansionRatio_;
```

Meaning: `expansionRatio` = end-size / start-size of the section. A negative value is treated as its inverse.

**File:** `gradingDescriptor/gradingDescriptor.C`

```cpp
Foam::gradingDescriptor::gradingDescriptor()
:
    blockFraction_(1),
    nDivFraction_(1),
    expansionRatio_(1)
{}

void Foam::gradingDescriptor::correct()
{
    if (expansionRatio_ < 0)
    {
        expansionRatio_ = 1.0/(-expansionRatio_);
    }
}

Foam::gradingDescriptor Foam::gradingDescriptor::inv() const
{
    return gradingDescriptor
    (
        blockFraction_,
        nDivFraction_,
        1.0/expansionRatio_
    );
}

Foam::Istream& Foam::operator>>(Istream& is, gradingDescriptor& gd)
{
    // Examine next token
    token t(is);

    if (t.isNumber())
    {
        gd.blockFraction_ = 1.0;
        gd.nDivFraction_ = 1.0;
        gd.expansionRatio_ = t.number();
    }
    else if (t.isPunctuation(token::BEGIN_LIST))
    {
        is >> gd.blockFraction_ >> gd.nDivFraction_ >> gd.expansionRatio_;
        is.readEnd("gradingDescriptor");
    }

    gd.correct();

    is.check(FUNCTION_NAME);
    return is;
}
```

### 4b. A list of sections: `gradingDescriptors`

**File:** `gradingDescriptor/gradingDescriptors.C`

```cpp
Foam::gradingDescriptors::gradingDescriptors()
:
    List<gradingDescriptor>(1, gradingDescriptor())
{}

void Foam::gradingDescriptors::normalise()
{
    scalar sumBlockFraction = 0;
    scalar sumNDivFraction = 0;

    for (const gradingDescriptor& gd : *this)
    {
        sumBlockFraction += gd.blockFraction_;
        sumNDivFraction += gd.nDivFraction_;
    }

    for (gradingDescriptor& gd : *this)
    {
        gd.blockFraction_ /= sumBlockFraction;
        gd.nDivFraction_  /= sumNDivFraction;
        gd.correct();
    }
}

Foam::gradingDescriptors Foam::gradingDescriptors::inv() const
{
    gradingDescriptors ret(*this);

    forAll(ret, i)
    {
        ret[i] = operator[](ret.size() - i - 1).inv();
    }

    return ret;
}

Foam::Istream& Foam::operator>>(Istream& is, gradingDescriptors& gds)
{
    // Examine next token
    token t(is);

    if (t.isNumber())
    {
        gds = gradingDescriptors(gradingDescriptor(t.number()));
        gds.correct();
    }
    else
    {
        is.putBack(t);

        // Read the list for gradingDescriptors
        is >> static_cast<List<gradingDescriptor>&>(gds);

        gds.normalise();
    }

    is.check(FUNCTION_NAME);
    return is;
}
```

Summary:

- `4` becomes a list with one section `{blockFraction 1, nDivFraction 1, ratio 4}`.
- `((0.2 0.3 4) (0.8 0.7 0.25))` becomes two sections, normalised so the fractions sum to 1.
- `inv()` at list level does two things: it reverses the section order and inverts each ratio.

---

## 5. Three values become twelve edges [SEEN]

**File:** `blockDescriptor/blockDescriptor.C`, `blockDescriptor::assignGradings`

```cpp
    switch (ratios.size())
    {
        case 0:
        {
            expand_.resize(12);
            expand_ = gradingDescriptors();
            break;
        }
        case 1:
        {
            // Identical in x/y/z-directions
            expand_.resize(12);
            expand_ = ratios[0];
            break;
        }
        case 3:
        {
            expand_.resize(12);

            // x-direction
            expand_[0]  = ratios[0];
            expand_[1]  = ratios[0];
            expand_[2]  = ratios[0];
            expand_[3]  = ratios[0];

            // y-direction
            expand_[4]  = ratios[1];
            expand_[5]  = ratios[1];
            expand_[6]  = ratios[1];
            expand_[7]  = ratios[1];

            // z-direction
            expand_[8]  = ratios[2];
            expand_[9]  = ratios[2];
            expand_[10] = ratios[2];
            expand_[11] = ratios[2];
            break;
        }
        case 12:
        {
            expand_ = ratios;
            break;
        }
        default:
        {
            ok = false;
            break;
        }
    }
```

For `(4 1 1)`: edges 0-3 (x) get `[4]`, edges 4-7 (y) get `[1]`, edges 8-11 (z) get `[1]`. So `simpleGrading` is `edgeGrading` with each value repeated 4 times.

---

## 6. Points are generated later: `createPoints()`

**File:** `blocks/block/blockCreate.C` [SEEN]

```cpp
void Foam::block::createPoints()
{
    // Set local variables for mesh specification
    const label ni = density().x();
    const label nj = density().y();
    const label nk = density().z();

    const point& p000 = blockPoint(0);
    const point& p100 = blockPoint(1);
    const point& p110 = blockPoint(2);
    const point& p010 = blockPoint(3);

    const point& p001 = blockPoint(4);
    const point& p101 = blockPoint(5);
    const point& p111 = blockPoint(6);
    const point& p011 = blockPoint(7);

    // List of edge point and weighting factors
    pointField p[12];
    scalarList w[12];
    const int nCurvedEdges = edgesPointsWeights(p, w);   // <-- calls step 7
```

Macros at the top of the file (so `w0` means "x-edge 0 at index i"):

```cpp
#define w0 w[0][i]
#define w1 w[1][i]
#define w2 w[2][i]
#define w3 w[3][i]

#define w4 w[4][j]
#define w5 w[5][j]
#define w6 w[6][j]
#define w7 w[7][j]

#define w8 w[8][k]
#define w9 w[9][k]
#define w10 w[10][k]
#define w11 w[11][k]
```

**[NOT SEEN] What triggers `block::createPoints()`.** `blockMesh::points()` calls `blockMesh::createPoints()` (in `blockMesh/blockMeshCreate.C`), and `block::createPoints()` is called from either the `block` constructor (`block.C`) or that function. I have not seen either file, so I cannot say whether `lineDivide` runs during construction or at the first `points()` call. To settle it:

```bash
grep -n "createPoints" $BM/blocks/block/*.C $BM/blocks/block/*.H $BM/blockMesh/*.C
```

---

## 7. Edge points and weights: where `lineDivide` is called [SEEN]

**File:** `blockDescriptor/blockDescriptorEdges.C`

### 7a. The loop over 12 edges

```cpp
int Foam::blockDescriptor::edgesPointsWeights
(
    pointField (&edgesPoints)[12],
    scalarList (&edgesWeights)[12]
) const
{
    int nCurved = 0;

    for (label edgei = 0; edgei < 12; ++edgei)  //< hexCell::nEdges()
    {
        nCurved += calcEdgePointsWeights
        (
            edgesPoints[edgei],
            edgesWeights[edgei],
            hexCell::modelEdges()[edgei],

            sizes()[edgei/4],   // 12 edges -> 3 components (x,y,z)
            expand_[edgei]
        );
    }

    return nCurved;
}
```

`sizes()[edgei/4]` gives N = 20 for x-edges (0-3), `ny` for 4-7, `nz` for 8-11.

### 7b. The three branches

```cpp
    for (const blockEdge& cedge : blockEdges_)
    {
        const int cmp = cedge.compare(thisEdge);

        if (cmp > 0)
        {
            // Curve has the same orientation

            // Divide the line
            const lineDivide divEdge(cedge, nDiv, expand);

            edgePoints  = divEdge.points();
            edgeWeights = divEdge.lambdaDivisions();

            return 1;  // Found curved-edge: done
        }
        else if (cmp < 0)
        {
            // Curve has the opposite orientation

            // Divide the line
            const lineDivide divEdge(cedge, nDiv, expand.inv());

            const pointField& p = divEdge.points();
            const scalarList& d = divEdge.lambdaDivisions();

            edgePoints.resize(p.size());
            edgeWeights.resize(d.size());

            // Copy in reverse order
            const label pn = (p.size() - 1);
            forAll(p, pi)
            {
                edgePoints[pi]  = p[pn - pi];
                edgeWeights[pi] = 1 - d[pn - pi];
            }

            return 1;  // Found curved-edge: done
        }
    }

    // Not curved-edge: divide the edge as a straight line

    // Get list of points for this block
    const pointField blockPoints(blockShape_.points(vertices_));

    lineDivide divEdge
    (
        blockEdges::lineEdge(blockPoints, cellModelEdge),
        nDiv,
        expand
    );

    edgePoints = divEdge.points();
    edgeWeights = divEdge.lambdaDivisions();

    return 0;
```

Your unit cube takes the last path (straight edge, `expand` used as is). `inv()` is only used for a curved edge stored in the opposite orientation.

---

## 8. THE FORMULA: `lineDivide` constructor [SEEN]

**File:** `blockEdges/lineDivide/lineDivide.C`

```cpp
namespace Foam
{
    //- Calculate the geometric expansion factor from the expansion ratio
    inline scalar calcGexp(const scalar expRatio, const label nDiv)
    {
        return nDiv > 1 ? pow(expRatio, 1.0/(nDiv - 1)) : 0.0;
    }
}

Foam::lineDivide::lineDivide
(
    const blockEdge& cedge,
    const label nDiv,
    const gradingDescriptors& gd
)
:
    points_(nDiv + 1),
    divisions_(nDiv + 1)
{
    divisions_[0]    = 0.0;
    divisions_[nDiv] = 1.0;

    scalar secStart = divisions_[0];
    label secnStart = 1;

    // Check that there are more divisions than sections
    if (nDiv >= gd.size())
    {
        // Calculate distribution of divisions to be independent
        // of the order of the sections
        labelList secnDivs(gd.size());
        label sumSecnDivs = 0;
        label secnMaxDivs = 0;

        forAll(gd, sectioni)
        {
            scalar nDivFrac = gd[sectioni].nDivFraction();
            secnDivs[sectioni] = label(nDivFrac*nDiv + 0.5);
            sumSecnDivs += secnDivs[sectioni];

            // Find the section with the largest number of divisions
            if (nDivFrac > gd[secnMaxDivs].nDivFraction())
            {
                secnMaxDivs = sectioni;
            }
        }

        // Adjust the number of divisions on the section with the largest
        // number so that the total is nDiv
        if (sumSecnDivs != nDiv)
        {
            secnDivs[secnMaxDivs] += (nDiv - sumSecnDivs);
        }

        forAll(gd, sectioni)
        {
            scalar blockFrac = gd[sectioni].blockFraction();
            scalar expRatio = gd[sectioni].expansionRatio();

            label secnDiv = secnDivs[sectioni];
            label secnEnd = secnStart + secnDiv;

            // Calculate the spacing
            if (equal(expRatio, 1))
            {
                for (label i = secnStart; i < secnEnd; i++)
                {
                    divisions_[i] =
                        secStart
                      + blockFrac*scalar(i - secnStart + 1)/secnDiv;
                }
            }
            else
            {
                // Calculate geometric expansion factor from the expansion ratio
                const scalar expFact = calcGexp(expRatio, secnDiv);

                for (label i = secnStart; i < secnEnd; i++)
                {
                    divisions_[i] =
                        secStart
                      + blockFrac*(1.0 - pow(expFact, i - secnStart + 1))
                    /(1.0 - pow(expFact, secnDiv));
                }
            }

            secStart = divisions_[secnEnd - 1];
            secnStart = secnEnd;
        }
    }
    // Otherwise mesh uniformly
    else
    {
        for (label i=1; i < nDiv; i++)
        {
            divisions_[i] = scalar(i)/nDiv;
        }
    }

    // Calculate the points
    points_ = cedge.position(divisions_);
}

const Foam::pointField& Foam::lineDivide::points() const noexcept
{
    return points_;
}

const Foam::scalarList& Foam::lineDivide::lambdaDivisions() const noexcept
{
    return divisions_;
}
```

### Derivation

Per-cell growth factor: `r = ratio^(1/(N-1))`, because the ratio is (last cell)/(first cell) = `r^(N-1)`.

Cell widths: `dx_i = dx_0 r^i`. Then `L = dx_0 (r^N - 1)/(r - 1)`, so `dx_0 = L (r - 1)/(r^N - 1)`.

Point locations (closed form, no running sum):

```
lambda_j = (1 - r^j) / (1 - r^N)        for one section (secStart = 0, blockFrac = 1)
```

Multiple sections: `secStart + blockFrac * (...)`, where each section starts at the end of the previous one.

### Worked numbers: N = 20, ratio 4

```
r        = 4^(1/19)        = 1.0757
r^20     = 4 * r           = 4.303
dx_0     = (r-1)/(r^20-1)  = 0.0757/3.303 = 0.0229     (uniform would be 0.05)
lambda_1 = 0.0229
dx_19    = dx_0 * r^19     = 4 * 0.0229    = 0.0917     (last/first = 4)
lambda_19 = 1 - 0.0917     = 0.9083
```

---

## 9. Lambda becomes a point on the edge

`points_ = cedge.position(divisions_);` (last line of the constructor). For a straight edge this is `lineEdge::position`, in `blockEdges/lineEdge/lineEdge.C` **[SEEN]**, see section 9b below.

```bash
sed -n '/position/,/^}/p' $BM/blockEdges/lineEdge/lineEdge.C
```

---

## 10. Interpolation into the 3D lattice [SEEN]

**File:** `blocks/block/blockCreate.C`, inner loop of `createPoints()`

```cpp
    points_.resize(nPoints());

    points_[pointLabel(0,  0,  0)] = p000;
    points_[pointLabel(ni, 0,  0)] = p100;
    points_[pointLabel(ni, nj, 0)] = p110;
    points_[pointLabel(0,  nj, 0)] = p010;
    points_[pointLabel(0,  0,  nk)] = p001;
    points_[pointLabel(ni, 0,  nk)] = p101;
    points_[pointLabel(ni, nj, nk)] = p111;
    points_[pointLabel(0,  nj, nk)] = p011;

    for (label k=0; k<=nk; k++)
    {
        for (label j=0; j<=nj; j++)
        {
            for (label i=0; i<=ni; i++)
            {
                // Skip block vertices
                if (vertex(i, j, k)) continue;

                const label vijk = pointLabel(i, j, k);

                // Calculate the weighting factors for all edges

                // x-direction
                scalar wx1 = (1 - w0)*(1 - w4)*(1 - w8) + w0*(1 - w5)*(1 - w9);
                scalar wx2 = (1 - w1)*w4*(1 - w11)      + w1*w5*(1 - w10);
                scalar wx3 = (1 - w2)*w7*w11            + w2*w6*w10;
                scalar wx4 = (1 - w3)*(1 - w7)*w8       + w3*(1 - w6)*w9;

                const scalar sumWx = wx1 + wx2 + wx3 + wx4;
                wx1 /= sumWx;
                wx2 /= sumWx;
                wx3 /= sumWx;
                wx4 /= sumWx;

                // (y- and z-direction weights wy1..wy4, wz1..wz4 are built the same way)

                // Points on straight edges
                const vector edgex1 = p000 + (p100 - p000)*w0;
                const vector edgex2 = p010 + (p110 - p010)*w1;
                const vector edgex3 = p011 + (p111 - p011)*w2;
                const vector edgex4 = p001 + (p101 - p001)*w3;
                // (edgey1..4 and edgez1..4 likewise, using w4..w7 and w8..w11)

                // Add the contributions
                points_[vijk] =
                (
                    wx1*edgex1 + wx2*edgex2 + wx3*edgex3 + wx4*edgex4
                  + wy1*edgey1 + wy2*edgey2 + wy3*edgey3 + wy4*edgey4
                  + wz1*edgez1 + wz2*edgez2 + wz3*edgez3 + wz4*edgez4
                )/3;

                // Apply curved-edge correction if block has curved edges
                if (nCurvedEdges)
                {
                    const vector corx1 = wx1*(p[0][i] - edgex1);
                    ...
                    points_[vijk] +=
                    (
                        corx1 + corx2 + corx3 + corx4
                      + cory1 + cory2 + cory3 + cory4
                      + corz1 + corz2 + corz3 + corz4
                    );
                }
            }
        }
    }

    if (!nCurvedFaces()) return;
    // ... three more loops apply face-curvature corrections (x-, y-, z-faces)
```

Key points:

- For straight edges, the mesh coordinates come from `w[e][i]` (the lambdas), not from `p[e][i]`. `p` is used only in the curved-edge correction. This is why changing `divisions_` in `lineDivide` is enough.
- Cell and boundary generation live in the same file:
    - `createCells()`: loops `k, j, i` (`i` fastest), `blockCells_[celli] = vertLabels(i, j, k)`.
    - `createBoundary()`: six patches in fixed order 0 = x-min, 1 = x-max, 2 = y-min, 3 = y-max, 4 = z-min, 5 = z-max.
    - `addBoundaryFaces()`: builds the 4 point labels of each boundary face.

---

## Appendix A: "modify" map

|Goal|File to copy and edit|
|---|---|
|Change the point distribution law|`blockEdges/lineDivide/lineDivide.C` (the `else` branch)|
|Change how the grading is written in the dict|`gradingDescriptor.C`, `gradingDescriptors.C`, `blockDescriptor.C` (`assignGradings`)|
|Change reversed-edge handling|`blockDescriptorEdges.C` (`cmp < 0` branch)|
|Change interior blending|`blocks/block/blockCreate.C` (`createPoints`)|
|Change edge shape|new class in `blockEdges/`|

Never edit `/usr/lib/openfoam`. Copy `src/mesh/blockMesh` to `~/OpenFOAM/dell-v2412/src/myBlockMesh`.

## Appendix B: `x = L*xi^p` patch for `lineDivide.C` (in your copy only)

Replace the geometric `else` branch. It is written so `inv()` (p -> 1/p) yields the mirror image, which the reversed-curved-edge branch in step 7b relies on.

```cpp
            else
            {
                // MODIFIED: power law. expRatio is now an exponent p.
                // p >= 1 : xi^p. p < 1 : mirror image, so inv() (p -> 1/p) works.
                for (label i = secnStart; i < secnEnd; i++)
                {
                    const scalar xi = scalar(i - secnStart + 1)/scalar(secnDiv);

                    const scalar f =
                        (expRatio >= 1)
                      ? pow(xi, expRatio)
                      : 1.0 - pow(1.0 - xi, 1.0/expRatio);

                    divisions_[i] = secStart + blockFrac*f;
                }
            }
```

Semantics changes to document:

- `4` no longer means last/first = 4. It means exponent 4. Use `(2 1 1)` for `x = L*xi^2`.
- `p = 1` is uniform, matching the existing `equal(expRatio, 1)` branch.
- A negative value is turned into `1/|p|` by `correct()`, giving the mirrored law.

Test table for `simpleGrading (2 1 1)`, N = 20, L = 1:

|j|expected x_j = (j/20)^2|
|---|---|
|1|0.0025|
|5|0.0625|
|10|0.25|
|19|0.9025|
|20|1.0|

Debug print (temporary, just before `points_ = cedge.position(divisions_);`):

```cpp
Info<< "lineDivide nDiv=" << nDiv << " lambda=" << divisions_ << endl;
```

## Appendix C: still to verify

|Item|Command|
|---|---|
|`block::iNew` forwarding|`grep -n "iNew" $BM/blocks/block/block.H $BM/blocks/block/block.C`|
|What calls `block::createPoints()`|`grep -n "createPoints" $BM/blocks/block/*.C $BM/blocks/block/*.H $BM/blockMesh/*.C`|
|`lineEdge::position`|`cat $BM/blockEdges/lineEdge/lineEdge.C`|
|Global assembly of block points|`cat $BM/blockMesh/blockMeshCreate.C`|
|Merging duplicate points|`cat $BM/blockMesh/blockMeshMergeTopological.C`|

---

## 9b. Lambda becomes a point on a straight edge [SEEN]

**File:** `blockEdges/lineEdge/lineEdge.C`

```cpp
namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(lineEdge, 0);
    addToRunTimeSelectionTable(blockEdge, lineEdge, Istream);
}
}

Foam::point Foam::blockEdges::lineEdge::position(const scalar lambda) const
{
    return blockEdge::linearPosition(lambda);
}

Foam::scalar Foam::blockEdges::lineEdge::length() const
{
    return Foam::mag(lastPoint() - firstPoint());
}
```

- `lineEdge::position` only delegates to `blockEdge::linearPosition(lambda)` (in `blockEdge.H`/`blockEdgeI.H`, **[NOT SEEN]**; expected to be `first + lambda*(last - first)`).
- `lineDivide` calls `position(const scalar)` through the `blockEdge` base type. For a curved edge, the same call goes to `arcEdge::position`, `splineEdge::position`, etc.
- `addToRunTimeSelectionTable(blockEdge, lineEdge, Istream)` is the registration pattern to copy for your own edge type.
- `lineDivide` is called with `blockEdges::lineEdge(blockPoints, cellModelEdge)`, the `(const pointField&, const edge&)` constructor, which is how step 7b builds a temporary straight edge.

---

## 11. Global assembly [SEEN]

**File:** `blockMesh/blockMeshCreate.C`

### 11a. `blockMesh::createPoints() const`: block points to global points

```cpp
    points_.resize(nPoints_);

    forAll(blocks, blocki)
    {
        const pointField& blockPoints = blocks[blocki].points();   // <-- already-computed block lattice

        const labelSubList pointAddr
        (
            mergeList_,
            blockPoints.size(),
            blockOffsets_[blocki]
        );

        UIndirectList<point>(points_, pointAddr) = blockPoints;
        ...
    }

    inplacePointTransforms(points_);   // scale / prescale / transform applied last
```

- `blocks[blocki].points()` is only read here. It is a `const` function, so the block's `points_` must already be filled. That strongly suggests `block::createPoints()` runs in the **block constructor** (during `createTopology`), not at the first `blockMesh::points()` call. This corrects what I said earlier. **Confirmed in section 13 (`block.C`).**
- `mergeList_` maps each block-local point (offset by `blockOffsets_[blocki]`) to its global point label, so duplicate points on shared block faces collapse into one. It is filled by `calcTopologicalMerge()` / `calcGeometricalMerge()`, in `blockMeshMergeTopological.C` / `blockMeshMergeGeometrical.C` **[NOT SEEN]**.
- `scale` is applied after assembly, so it never changes grading, only units.

### 11b. Verbose cell-size printout: a free test for grading

With `verbose_` on (the default), `createPoints` prints for each block the first and last cell width along i, j, k, measured on the lines through vertex `(0,0,0)`, multiplied by the scale:

```cpp
const scalar cwBeg = mag(blockPoints[v1] - blockPoints[v0]);
const scalar cwEnd = mag(blockPoints[vn] - blockPoints[vn1]);
Info<< "        i : " << cwBeg*scaleTot.x() << " .. " << cwEnd*scaleTot.x() << nl;
```

Expected output for `(20 20 20) simpleGrading (4 1 1)` on a unit cube:

```
    Block 0 cell size :
        i : 0.0229... .. 0.0917...
        j : 0.05 .. 0.05
        k : 0.05 .. 0.05
```

Expected for `simpleGrading (2 1 1)` with the power-law patch: `i : 0.0025 .. 0.0975`, since `(1/20)^2 = 0.0025` and `1 - (19/20)^2 = 0.0975`.

### 11c. Cells, patch faces and the final polyMesh

```cpp
void Foam::blockMesh::createCells() const
{
    ...
    forAll(blocks, blocki)
    {
        for (const hexCell& blockCell : blocks[blocki].cells())
        {
            forAll(cellPoints, cellPointi)
            {
                cellPoints[cellPointi] =
                    mergeList_[ blockCell[cellPointi] + blockOffsets_[blocki] ];
            }

            // Construct collapsed cell and add to list
            cells_[celli].reset(hex, cellPoints, true);
            ++celli;
        }
    }
}
```

- Cells: each block cell's 8 local point labels are mapped through `mergeList_` to global labels. The `true` in `reset(..., true)` allows collapsed (degenerate) hex cells.
- `createPatchFaces(polyPatch)`: for each topology patch face, finds the matching block face, maps the boundary quads through `mergeList_`, and **collapses duplicate point labels**. A quad with 4 unique points is kept as a quad, one with 3 becomes a triangle, and a face collapsed to an edge or point is dropped.
- `createPatches()`: calls `createPatchFaces` for every topology patch.
- `blockMesh::mesh(io)`: builds the final `polyMesh` from `points()`, `cells()`, `patches()`, `patchNames()`, `patchDicts()`, with default patch `"defaultFaces"` of type `empty`. It also creates cell zones from named blocks. Merge patch pairs and cyclics are done elsewhere (the application's `.H` snippets, which need `libdynamicMesh`).

---

## Updated pipeline (all steps now seen except the ones marked)

```
 1  blockMesh ctor -> createTopology()                         blockMesh.C
 2  createTopology(): blocks via block::iNew                   blockMeshTopology.C
      -> block ctor -> (probably) block::createPoints()        block.C [SEEN, section 13]
 3  blockDescriptor ctor parses hex, N, List<gradingDescriptors>  blockDescriptor.C
 4  gradingDescriptors operator>>: 4 -> [{1,1,4}]              gradingDescriptors.C
 5  assignGradings(): 12 edges                                 blockDescriptor.C
 6  createPoints(): edgesPointsWeights(p, w)                   blockCreate.C
 7  calcEdgePointsWeights(): lineDivide(edge, N, expand)       blockDescriptorEdges.C
 8  lineDivide ctor: lambda_j formula                          lineDivide.C
 9  points_ = cedge.position(lambda) -> linearPosition         lineEdge.C, blockEdge.C [SEEN, section 13]; linearPosition in blockEdgeI.H [NOT SEEN]
10  w[e][i]=lambda -> blended -> block points_                 blockCreate.C
11  blockMesh::createPoints(): block points -> global via mergeList_, then scale
                                                               blockMeshCreate.C
12  blockMesh::mesh(): polyMesh(points, cells, patches)        blockMeshCreate.C
```

Still open: `block.H` (`block::iNew`), `blockEdgeI.H` (`linearPosition`).

---

## 12. Merging duplicate points: `calcTopologicalMerge()` [SEEN]

**File:** `blockMesh/blockMeshMergeTopological.C`, called from the `blockMesh` constructor (default strategy `MERGE_TOPOLOGY`). It runs **after** all blocks, and therefore all `lineDivide` calls, are done.

### 12a. Offsets: every block gets a range of global point slots

```cpp
    forAll(blocks, blocki)
    {
        blockOffsets_[blocki] = nPoints_;

        nPoints_ += blocks[blocki].nPoints();
        nCells_  += blocks[blocki].nCells();
    }
    mergeList_.setSize(nPoints_, -1);
```

Block-local point `p` of block `b` has "unmerged" global index `p + blockOffsets_[b]`.

### 12b. Match points on each shared (internal) block face

Only the internal faces of the _topology_ mesh (the coarse mesh whose cells are the blocks) are processed. For each one, the owner block P and neighbour block N are found, and `faceMap(...)` works out how the (i, j) axes of the two faces correspond (rotation and flip), using the static tables `faceEdgeDirs` and `faceFaceRotMap`.

The grading consistency checks are here:

```cpp
        // Check block subdivision correspondence
        ...
            if (Pnij != NPnij)
            {
                FatalErrorInFunction
                    << "Sub-division mismatch between face " ...
            }

        const boundBox bb(topoCells[blockPi].points(topoFaces, topoPoints));
        const scalar testSqrDist = magSqr(1e-6*bb.span());
        ...
                if (sqrDist > testSqrDist)
                {
                    FatalErrorInFunction
                        << "Point merge failure between face " ...
                        << "    This may be due to inconsistent grading."
                        << exit(FatalError);
                }
```

- Cell counts on the shared face must match, or you get "Sub-division mismatch".
- Every matched pair of points must lie within `1e-6 * block span` of each other, or you get "Point merge failure ... inconsistent grading".

**What this means for your `x = L*xi^p` law:** two neighbouring blocks that share an x-edge, but traverse it in opposite directions, must give the _same physical_ points. With geometric grading, that is why the neighbour uses ratio `1/r`. The mirror-aware patch (Appendix B) keeps this property: block P with `p` gives `xi^p`, and block N traversing backwards with `1/p` gives `1 - (1-xi)^p`, which is the same point set. A plain `xi^p` law without the mirror branch would fail this check with "inconsistent grading".

### 12c. Union of matched points (smallest label wins)

```cpp
                label Ppointi = blockPpointi + blockOffsets_[blockPi];
                label Npointi = blockNpointi + blockOffsets_[blockNi];

                label minPNi = min(Ppointi, Npointi);

                if (mergeList_[Ppointi] != -1)
                {
                    minPNi = min(minPNi, mergeList_[Ppointi]);
                }

                if (mergeList_[Npointi] != -1)
                {
                    minPNi = min(minPNi, mergeList_[Npointi]);
                }

                mergeList_[Ppointi] = mergeList_[Npointi] = minPNi;
```

Then a `do { ... } while (changedPointMerge)` loop repeats the same pass over all faces until nothing changes (needed for points on shared edges and corners that belong to more than two blocks). It aborts after 100 passes.

### 12d. Compaction into consecutive global labels

```cpp
    label nUniqPoints = 0;

    forAll(mergeList_, pointi)
    {
        if (mergeList_[pointi] == -1 || mergeList_[pointi] == pointi)
        {
            mergeList_[pointi] = nUniqPoints++;      // a new unique point
        }
        else
        {
            mergeList_[pointi] = mergeList_[mergeList_[pointi]];   // copy the label of the point it merged into
        }
    }

    nPoints_ = nUniqPoints;
```

After this, `mergeList_[unmerged index]` is the final global point label used by `blockMesh::createPoints()`, `createCells()` and `createPatchFaces()` (section 11).

`MERGE_POINTS` (`calcGeometricalMerge`, in `blockMeshMergeGeometrical.C`, **[NOT SEEN]**) fills the same `mergeList_` by distance instead. It is chosen automatically if `checkDegenerate()` finds collapsed blocks, or via `mergeType points;`.

---

## 13. `block.C` and `blockEdge.C` [SEEN]

### 13a. Answer to "when does `lineDivide` run": in the `block` constructor

**File:** `blocks/block/block.C`

```cpp
Foam::block::block
(
    const dictionary& dict,
    const label index,
    const pointField& vertices,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    Istream& is
)
:
    blockDescriptor(dict, index, vertices, edges, faces, is),   // parses (20 20 20) simpleGrading (4 1 1)
    points_(),
    blockCells_(),
    blockPatches_()
{
    // Always need points, and demand-driven data leaves dangling addressing?
    createPoints();      // -> edgesPointsWeights -> lineDivide (steps 6-10)
    createBoundary();
}
```

The same two calls are in the other two constructors. So the corrected timeline is:

```
blockMesh ctor
  -> createTopology()
       -> block::iNew / block::New  -> block ctor
            -> blockDescriptor ctor : parse grading, assignGradings   (steps 3-5)
            -> createPoints()       : lineDivide + interpolation      (steps 6-10)
            -> createBoundary()
  -> calcTopologicalMerge()        (section 12, uses the block points)
later: points() -> blockMesh::createPoints() copies block points to global (section 11)
```

`createCells()` for a block is not called here; block cells are created on demand (`blockCells_` filled when needed). Only the global `blockMesh::points()`, `cells()` and `patches()` are lazy.

### 13b. Block selector

```cpp
Foam::autoPtr<Foam::block> Foam::block::New(dict, index, points, edges, faces, Istream& is)
{
    const word blockOrCellShapeType(is);

    auto* ctorPtr = IstreamConstructorTable(blockOrCellShapeType);

    if (!ctorPtr)
    {
        is.putBack(token(blockOrCellShapeType));
        return autoPtr<block>::New(dict, index, points, edges, faces, is);
    }

    return autoPtr<block>(ctorPtr(dict, index, points, edges, faces, is));
}
```

It reads the first word (`hex`). If it is not a registered block type, the word is put back and a normal `block` is built. So `hex` is not a registered type, and the default `block` constructor above reads it as the cell model. The `iNew` wrapper in `block.H` **[NOT SEEN]** presumably calls this `New`.

### 13c. `blockEdge`: base class of every edge

**File:** `blockEdges/blockEdge/blockEdge.C`

```cpp
Foam::blockEdge::blockEdge(const pointField& points, const edge& fromTo)
:
    points_(points),
    start_(fromTo.first()),
    end_(fromTo.last())
{}

Foam::tmp<Foam::pointField>
Foam::blockEdge::position(const scalarList& lambdas) const
{
    auto tpoints = tmp<pointField>::New(lambdas.size());
    auto& points = tpoints.ref();

    forAll(lambdas, i)
    {
        points[i] = position(lambdas[i]);
    }
    return tpoints;
}
```

- **This is the exact call in `lineDivide`:** `points_ = cedge.position(divisions_);` uses the `scalarList` overload, which loops over the lambdas and calls the virtual `position(const scalar)` of the concrete edge (`lineEdge`, `arcEdge`, ...). For `lineEdge` that is `blockEdge::linearPosition(lambda)` (in `blockEdgeI.H`, **[NOT SEEN]**).
- The edge stores only vertex _indices_ (`start_`, `end_`) and a reference to the point array `points_`.
- The runtime selector: `blockEdge::New` reads the type word (`arc`, `spline`, `line`, ...) and looks it up in `IstreamConstructorTable`. An unknown word gives `FatalIOErrorInLookup`. This is what your custom edge library plugs into.

### 13d. Where `compare()` comes from

`calcEdgePointsWeights` calls `cedge.compare(thisEdge)`. It is declared in the `blockEdge` header/inline file **[NOT SEEN]**. It returns `> 0` for the same orientation, `< 0` for opposite, `0` for no match.