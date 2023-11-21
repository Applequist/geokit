use std::fmt;

/// A pair that uniquely identified a CRS element.
pub type AuthorityId = (String, u16);

/// An CRS element id.
#[derive(Debug, Clone)]
pub struct Id {
    name: String,
    id: Option<AuthorityId>,
}

impl Id {
    pub fn name<S: Into<String>>(name: S) -> Self {
        Self {
            name: name.into(),
            id: None,
        }
    }

    pub fn full<S: Into<String>>(name: S, authority: S, code: u16) -> Self {
        Self {
            name: name.into(),
            id: Some((authority.into(), code)),
        }
    }

    pub fn renamed<S: Into<String>>(&self, rename: S) -> Self {
        Self {
            name: rename.into(),
            id: self.id.clone(),
        }
    }
}

impl From<&str> for Id {
    fn from(value: &str) -> Self {
        Id::name(value)
    }
}

impl From<(&str, u16)> for Id {
    fn from(value: (&str, u16)) -> Self {
        Id::full("n/a", value.0, value.1)
    }
}

impl PartialEq for Id {
    fn eq(&self, other: &Self) -> bool {
        match (&self.id, &other.id) {
            (Some((this_auth, this_code)), Some((other_auth, other_code))) => {
                this_auth == other_auth && this_code == other_code
            }
            (Some(_), None) | (None, Some(_)) => false,
            _ => self.name == other.name,
        }
    }
}

impl Eq for Id {}

impl fmt::Display for Id {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(authority) = &self.id {
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
        // If no id, consider name only:
        let name = Id::name("WGS 84");
        assert!(name.eq(&name));
        assert!(!name.ne(&name));
        assert!(name.eq(&"WGS 84".into()));
        assert!(!name.ne(&"WGS 84".into()));

        let other_name = Id::name("WGS 84.1");
        assert!(!name.eq(&other_name));
        assert!(!other_name.eq(&name));
        assert!(name.ne(&other_name));
        assert!(other_name.ne(&name));

        // If id present, only consider id
        let id = Id::full("WGS 84", "EPSG", 7030);
        assert!(id.eq(&id));
        assert!(!id.ne(&id));
        assert!(!name.eq(&id));
        assert!(name.ne(&id));
        assert!(!id.eq(&name));
        assert!(id.ne(&name));

        let same_id = Id::full("WGS 1984", "EPSG", 7030);
        assert!(id.eq(&same_id));
        assert!(!id.ne(&same_id));
        assert!(same_id.eq(&id));
        assert!(!same_id.ne(&id));

        let other_id = Id::full("WGS 84", "geokit", 7030);
        assert!(!id.eq(&other_id));
        assert!(id.ne(&other_id));
        assert!(!other_id.eq(&id));
        assert!(other_id.ne(&id));
    }

    #[test]
    fn display() {
        let id = Id::name("WGS 84");
        assert_eq!(id.to_string(), String::from("WGS 84"));

        let full = Id::full("WGS 84", "epsg", 7030);
        assert_eq!(full.to_string(), "WGS 84 (epsg:7030)");
    }
}
