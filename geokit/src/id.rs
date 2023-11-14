use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Id {
    name: String,
    authority: Option<(String, u16)>,
}

impl Id {
    pub fn name<S: Into<String>>(name: S) -> Self {
        Self {
            name: name.into(),
            authority: None,
        }
    }

    pub fn full<S: Into<String>>(name: S, authority: S, code: u16) -> Self {
        Self {
            name: name.into(),
            authority: Some((authority.into(), code)),
        }
    }
}

impl fmt::Display for Id {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(authority) = &self.authority {
            write!(f, "{} ({}:{})", self.name, authority.0, authority.1)
        } else {
            write!(f, "{}", self.name)
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn name() {
        let id = Id::name("WGS 84");
        assert_eq!(id.to_string(), String::from("WGS 84"));
    }

    #[test]
    fn full() {
        let full_id = Id::full("WGS 84", "epsg", 7030);
        assert_eq!(full_id.to_string(), "WGS 84 (epsg:7030)");
    }

    #[test]
    fn clone() {
        let id = Id::name("WGS 84");
        let cpy = id.clone();
        assert_eq!(id, cpy);

        let full = Id::full("WGS 84", "espg", 7030);
        let full_cpy = full.clone();
        assert_eq!(full, full_cpy);
    }

    #[test]
    fn eq() {
        let id = Id::name("WGS 84");
        assert!(id.eq(&id));
        assert!(!id.ne(&id));

        let id2 = Id::name("WGS 84.1");
        assert!(!id.eq(&id2));
        assert!(!id2.eq(&id));
        assert!(id.ne(&id2));
        assert!(id2.ne(&id));

        let full = Id::full("WGS 84", "epsg", 7030);
        assert!(!id.eq(&full));
        assert!(id.ne(&full));
        assert!(!full.eq(&id));
        assert!(full.ne(&id));

        let full2 = Id::full("WGS 84", "geokit", 7030);
        assert!(!full.eq(&full2));
        assert!(full.ne(&full2));
        assert!(!full2.eq(&full));
        assert!(full2.ne(&full));
    }

    #[test]
    fn display() {
        let id = Id::name("WGS 84");
        assert_eq!(id.to_string(), String::from("WGS 84"));

        let full = Id::full("WGS 84", "epsg", 7030);
        assert_eq!(full.to_string(), "WGS 84 (epsg:7030)");
    }
}
