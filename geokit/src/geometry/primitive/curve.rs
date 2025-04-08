use super::{Boundary, point::Point};
use crate::{
    crs::Crs,
    geometry::{
        Geometry, GeometryType,
        complex::Complex,
        coordinate::{
            Pos,
            curve::{CurveSegment, ParameterizedCurve, line_string::LineString},
        },
    },
    quantities::length::Length,
    units::length::M,
};
use dyn_clone::clone_box;
use std::fmt::Debug;

/// A [Boundary] used to describe a [Curve](crate::geometry::primitive::curve)'s boundary if any.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CurveBoundary {
    pub start: Point,
    pub end: Point,
}

impl Geometry for CurveBoundary {
    fn geometry_type(&self) -> GeometryType {
        self.start.geometry_type()
    }

    fn is_cycle(&self) -> bool {
        true
    }

    fn boundary(&self) -> Option<Box<dyn Boundary>> {
        None
    }
}

impl Complex for CurveBoundary {}

impl Boundary for CurveBoundary {}

/// [Curve] is the 1-dimensional primitive geometric object.
///
/// A [Curve] is the image of an **open** interval into a 2- or 3-dimensional Euclidian space
/// by a continuous mapping. As such it can be defined as a [parameterized
/// function](crate::geometry::coordinate::curve::ParameterizedCurve)
///
/// The internal representation of a curve is a sequence of [CurveSegment] with the appropriate
/// start/end relationship.
///
/// A [Curve] is created via a [CurveBuilder] to ensure these relationships are satisfied.
#[derive(Clone, Debug)]
pub struct Curve {
    coord_dim: usize,
    seg: Vec<Box<dyn CurveSegment>>,
}

impl Geometry for Curve {
    fn geometry_type(&self) -> GeometryType {
        GeometryType::Curve(self.coord_dim)
    }

    fn is_cycle(&self) -> bool {
        self.start() == self.end()
    }

    fn boundary(&self) -> Option<Box<dyn Boundary>> {
        if self.is_cycle() {
            None
        } else {
            Some(Box::new(CurveBoundary {
                start: Point::new(self.start()),
                end: Point::new(self.end()),
            }))
        }
    }
}

impl ParameterizedCurve for Curve {
    fn coord_dim(&self) -> usize {
        self.coord_dim
    }

    fn start(&self) -> &Pos {
        self.seg.first().unwrap().start()
    }

    fn end(&self) -> &Pos {
        self.seg.last().unwrap().end()
    }

    fn length(&self, crs: &dyn Crs) -> Length {
        self.seg
            .iter()
            .map(|s| s.length(crs))
            .fold(0. * M, |acc, l| acc + l)
    }

    fn param(&self, crs: &dyn Crs, s: Length) -> Box<Pos> {
        assert!(s <= self.length(crs), "Invalid curvilinear parameter `s`");
        let mut seg_iter = self.seg.iter();
        let mut res_len = s;
        let mut seg_ix = 0;
        while let Some(s) = seg_iter.next() {
            let seg_len = s.length(crs);
            if seg_len >= res_len {
                break;
            }
            res_len -= seg_len;
            seg_ix += 1;
        }
        let seg = &self.seg[seg_ix];
        seg.param(crs, res_len)
    }

    fn as_line_string(
        &self,
        crs: &dyn Crs,
        max_distance: Option<Length>,
        max_offset: Option<Length>,
    ) -> LineString {
        todo!()
    }
}

#[derive(Clone, Debug)]
pub struct CurveBuilder<'a> {
    crs: &'a dyn Crs,
    coord_dim: usize,
    seg: Vec<Box<dyn CurveSegment>>,
}

impl<'a> CurveBuilder<'a> {
    pub(crate) fn new(crs: &'a dyn Crs, seg: Vec<Box<dyn CurveSegment>>) -> Self {
        Self {
            crs,
            coord_dim: crs.dim(),
            seg,
        }
    }

    pub fn with_seg<C: CurveSegment + 'static>(
        crs: &'a dyn Crs,
        seg: C,
    ) -> Result<Self, &'static str> {
        Ok(Self::new(crs, vec![clone_box(&seg)]))
    }

    pub fn append_segment<C: CurveSegment + 'static>(
        mut self,
        seg: C,
    ) -> Result<Self, &'static str> {
        if self.coord_dim != seg.coord_dim() {
            Err("Invalid coordinate dimension")
        } else if self.seg.last().unwrap().end() != seg.start() {
            Err("Disconnected segment")
        } else {
            self.seg.push(clone_box(&seg));
            Ok(self)
        }
    }

    pub fn build(self) -> Curve {
        Curve {
            coord_dim: self.coord_dim,
            seg: self.seg,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        crs::projected::ProjectedCrs,
        cs::cartesian::projected::ProjectedAxes,
        geodesy::geodetic_datum::consts::WGS84,
        geometry::{
            Geometry, GeometryType,
            coordinate::curve::{ParameterizedCurve, line_string::LineStringBuilder},
        },
        projections::ProjectionSpec,
        units::length::M,
    };

    use super::CurveBuilder;

    #[test]
    fn curve() -> Result<(), &'static str> {
        let crs = ProjectedCrs {
            id: "UTM Zone 1".into(),
            datum: WGS84,
            axes: ProjectedAxes::EastNorth { horiz_unit: M },
            projection: ProjectionSpec::UTMNorth { zone: 1 },
        };

        let curve = CurveBuilder::with_seg(
            &crs,
            LineStringBuilder::with_line(&crs, &[0., 0.], &[1., 0.])?
                .line_to([1., 1.])?
                .build(),
        )?
        .append_segment(
            LineStringBuilder::with_line(&crs, &[1., 1.], &[0., 1.])?
                .line_to([0., 0.])?
                .build(),
        )?
        .build();

        assert_eq!(curve.geometry_type(), GeometryType::Curve(2));

        //assert!(
        //    (CurveBoundary {
        //        start: Point::new(curve.start()),
        //        end: Point::new(curve.end()),
        //    } as Boundary)
        //        .eq(curve.boundary().as_ref())
        //);
        assert!(curve.is_cycle());
        assert_eq!(curve.length(&crs), 4. * M);

        Ok(())
    }
}
