use super::{CurveInterpolation, CurveSegment, ParameterizedCurve};
use crate::{
    crs::Crs,
    geometry::{
        coordinate::Pos,
        primitive::{curve::CurveBoundary, point::Point},
    },
    math::fp::Float,
    quantities::length::Length,
    units::length::M,
};

/// [LineString] is a linear [CurveSegment].
#[derive(Clone, Debug, PartialEq)]
pub struct LineString {
    coord_dim: usize,
    pos: Vec<Float>,
    len: Vec<Length>,
}

impl LineString {
    /// Creates new line string.
    ///
    /// The following must hold true:
    /// - `pos.len() % coord_dim == 0`, i.e all control point positions must be *complete*.
    /// - `pos.len() >= 2 * coord_dim`, i.e at least 2 control points must be given.
    pub(crate) fn new(coord_dim: usize, pos: Vec<Float>, len: Vec<Length>) -> Self {
        assert!(pos.len() % coord_dim == 0, "Incomplete position in `pos`");
        assert!(
            pos.len() >= 2 * coord_dim,
            "Need at least 2 positions in `pos`"
        );
        assert!(
            len.len() == pos.len() / coord_dim - 1,
            "Invalid lengthes in `len`"
        );
        Self {
            coord_dim,
            pos,
            len,
        }
    }

    /// Returns the number of control position of this line string.
    pub(crate) fn len(&self) -> usize {
        self.pos.len() / self.coord_dim
    }

    fn pos_start(&self, n: usize) -> usize {
        n * self.coord_dim
    }

    /// Returns the coordinates of the `n`-th control position of this line string.
    pub(crate) fn pos(&self, n: usize) -> &Pos {
        &self.pos[self.pos_start(n)..self.pos_start(n + 1)]
    }

    pub(crate) fn pos_iter(&self) -> impl Iterator<Item = &Pos> {
        PosIter {
            ls: self,
            current: 0,
        }
    }
}

struct PosIter<'a> {
    ls: &'a LineString,
    current: usize,
}

impl<'a> Iterator for PosIter<'a> {
    type Item = &'a Pos;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.ls.len() {
            self.current += 1;
            Some(self.ls.pos(self.current - 1))
        } else {
            None
        }
    }
}

impl ParameterizedCurve for LineString {
    fn coord_dim(&self) -> usize {
        self.coord_dim
    }

    fn start(&self) -> &Pos {
        self.pos(0)
    }

    fn end(&self) -> &Pos {
        self.pos(self.len() - 1)
    }

    fn length(&self, _crs: &dyn Crs) -> Length {
        self.len.iter().fold(0. * M, |acc, l| acc + *l)
    }

    fn param(&self, _crs: &dyn Crs, s: Length) -> Box<Pos> {
        assert!(s <= self.length(_crs), "Invalid curvilinear parameter `s`");
        let mut len_iter = self.len.iter();
        let mut res_len = s;
        let mut start_pos_ix = 0;
        let mut ratio = 0.;
        while let Some(l) = len_iter.next() {
            if l >= &res_len {
                ratio = res_len / *l;
                break;
            }
            res_len -= *l;
            start_pos_ix += 1;
        }
        let start = self.pos(start_pos_ix);
        let end = self.pos(start_pos_ix + 1);

        start
            .iter()
            .zip(end.iter())
            .map(|(&s, &e)| (1. - ratio) * s + ratio * e)
            .collect::<Vec<_>>()
            .into_boxed_slice()
    }

    fn as_line_string(
        &self,
        crs: &dyn Crs,
        max_distance: Option<Length>,
        max_offset: Option<Length>,
    ) -> LineString {
        // TODO: don't forget to take max_distance in consideration!!
        todo!()
    }
}

impl CurveSegment for LineString {
    fn interpolation(&self) -> CurveInterpolation {
        CurveInterpolation::Linear
    }
    fn curve_boundary(&self) -> CurveBoundary {
        CurveBoundary {
            start: Point::new(self.pos(0)),
            end: Point::new(self.pos(self.len() - 1)),
        }
    }
}

/// [LineStringBuilder] helps creating [LineString].
#[derive(Clone, Debug)]
pub struct LineStringBuilder<'a> {
    crs: &'a dyn Crs,
    coord_dim: usize,
    pos: Vec<Float>,
    len: Vec<Length>,
}

impl<'a> LineStringBuilder<'a> {
    /// Creates a new [LineStringBuilder] with the given [Crs] and positions.
    ///
    /// The following must hold true:
    /// - `pos.len() % crs.dim() == 0`, i.e. all positions are completely specified,
    /// - `pos.len() >= 2 * crs.dim()`, i.e. at least 2 positions are given.
    /// - `len.len() == pos.len() / crs.dim() - 1, i.e.
    pub(crate) fn new(
        crs: &'a dyn Crs,
        pos: Vec<Float>,
        len: Vec<Length>,
    ) -> Result<Self, &'static str> {
        let coord_dim = crs.dim();
        assert!(pos.len() % coord_dim == 0, "Incomplete position in `pos`");
        assert!(
            pos.len() >= 2 * coord_dim,
            "Need at least 2 positions in `pos`"
        );
        assert!(
            len.len() == pos.len() / coord_dim - 1,
            "Invalid lengthes in `len`"
        );
        Ok(Self {
            crs,
            coord_dim,
            pos,
            len,
        })
    }

    /// Creates a new [LineStringBuilder] with the initial line
    /// segment from `start` to `to`.
    pub fn with_line<P: AsRef<Pos>>(
        crs: &'a dyn Crs,
        from: P,
        to: P,
    ) -> Result<Self, &'static str> {
        let mut pos = vec![];
        let from_pos = from.as_ref();
        let to_pos = to.as_ref();
        pos.extend(from_pos);
        pos.extend(to_pos);
        let len = vec![crs.dist(from_pos, to_pos)?];
        Self::new(crs, pos, len)
    }

    fn pos_start(&self, n: usize) -> usize {
        n * self.coord_dim
    }

    fn pos(&self, n: usize) -> &Pos {
        &self.pos[self.pos_start(n)..self.pos_start(n + 1)]
    }

    fn len(&self) -> usize {
        self.pos.len() / self.coord_dim
    }

    /// Adds a new line segment from the end point of the last segment
    /// to `p`.
    ///
    /// `p.as_ref().len() == self.coord_dim`
    pub fn line_to<P: AsRef<Pos>>(mut self, p: P) -> Result<Self, &'static str> {
        let s = p.as_ref();
        let to = (s.len() == self.coord_dim)
            .then_some(s)
            .ok_or("Invalid coordinate dimension")?;

        let length = self.crs.dist(self.pos(self.len() - 1), to)?;
        self.pos.extend(s);
        self.len.push(length);
        Ok(self)
    }

    /// Adds a copy of the start position to close the [LineString]
    /// and make it.
    pub fn close(self) -> Result<LineString, &'static str> {
        let start = self.pos(0).to_vec();
        self.line_to(start).map(|builder| builder.build())
    }

    /// Returns the [LineString] made from this builder.
    pub fn build(self) -> LineString {
        LineString {
            coord_dim: self.coord_dim,
            pos: self.pos,
            len: self.len,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        crs::projected::ProjectedCrs, cs::cartesian::projected::ProjectedAxes,
        geodesy::geodetic_datum::consts::WGS84, geometry::coordinate::curve::ParameterizedCurve,
        projections::ProjectionSpec, units::length::M,
    };

    use super::LineStringBuilder;

    #[test]
    fn line_string_pos() -> Result<(), &'static str> {
        let crs = ProjectedCrs {
            id: "UTM Zone 1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorth { horiz_unit: M },
            projection: ProjectionSpec::UTMNorth { zone: 1 },
        };

        let line_string = LineStringBuilder::with_line(&crs, &[0., 0.], &[1., 0.])?
            .line_to([1., 1.])?
            .line_to(vec![0., 1.])?
            .close()?;

        assert_eq!(line_string.len(), 5);
        assert_eq!(line_string.coord_dim, 2);
        assert_eq!(line_string.pos(0), [0., 0.]);
        assert_eq!(line_string.pos(2), [1., 1.]);
        assert_eq!(line_string.pos(4), [0., 0.]);

        Ok(())
    }

    #[test]
    fn line_string_pos_iter() -> Result<(), &'static str> {
        let crs = ProjectedCrs {
            id: "UTM Zone 1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorth { horiz_unit: M },
            projection: ProjectionSpec::UTMNorth { zone: 1 },
        };

        let line_string = LineStringBuilder::with_line(&crs, &[0., 0.], &[1., 0.])?
            .line_to([1., 1.])?
            .line_to(vec![0., 1.])?
            .close()?;

        let mut pos_iter = line_string.pos_iter();
        assert_eq!(pos_iter.next(), Some(&[0., 0.][..]));
        pos_iter.next();
        assert_eq!(pos_iter.next(), Some(&[1., 1.][..]));
        pos_iter.next();
        assert_eq!(pos_iter.next(), Some(&[0., 0.][..]));
        assert_eq!(pos_iter.next(), None);

        Ok(())
    }

    #[test]
    fn parameterized_curve() -> Result<(), &'static str> {
        let crs = ProjectedCrs {
            id: "UTM 1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorth { horiz_unit: M },
            projection: ProjectionSpec::UTMNorth { zone: 1 },
        };

        let line_string = LineStringBuilder::with_line(&crs, &[0., 0.], &[1., 0.])?
            .line_to([1., 1.])?
            .line_to(vec![0., 1.])?
            .close()?;

        assert_eq!(line_string.coord_dim(), 2);
        assert_eq!(line_string.start(), [0., 0.]);
        assert_eq!(line_string.end(), [0., 0.]);
        assert_eq!(line_string.length(&crs), 4. * M);
        assert_eq!(
            line_string.param(&crs, 0.5 * M),
            vec![0.5, 0.].into_boxed_slice()
        );

        Ok(())
    }
}
