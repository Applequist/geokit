use std::fmt;

/// A pair that uniquely identified a CRS element.
pub type AuthorityId = (String, u16);

/// An [Tag] adds a mandatory name and an optional id to an element.
#[derive(Debug, Clone)]
pub struct Tag {
    name: String,
    id: Option<AuthorityId>,
}

impl Tag {
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

impl From<&str> for Tag {
    fn from(value: &str) -> Self {
        Tag::name(value)
    }
}

impl From<(&str, &str, u16)> for Tag {
    fn from(value: (&str, &str, u16)) -> Self {
        Tag::full(value.0, value.1, value.2)
    }
}

impl PartialEq for Tag {
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

impl Eq for Tag {}

impl fmt::Display for Tag {
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
        let tag = Tag::name("WGS 84");
        assert_eq!(tag.to_string(), String::from("WGS 84"));
    }

    #[test]
    fn full() {
        let full_tag = Tag::full("WGS 84", "epsg", 7030);
        assert_eq!(full_tag.to_string(), "WGS 84 (epsg:7030)");
    }

    #[test]
    fn clone() {
        let tag = Tag::name("WGS 84");
        let cpy = tag.clone();
        assert_eq!(tag, cpy);

        let full_tag = Tag::full("WGS 84", "espg", 7030);
        let full_cpy = full_tag.clone();
        assert_eq!(full_tag, full_cpy);
    }

    #[test]
    fn eq() {
        // If no id, consider name only:
        let name_tag = Tag::name("WGS 84");
        assert!(name_tag.eq(&name_tag));
        assert!(!name_tag.ne(&name_tag));
        assert!(name_tag.eq(&"WGS 84".into()));
        assert!(!name_tag.ne(&"WGS 84".into()));

        let other_name_trag = Tag::name("WGS 84.1");
        assert!(!name_tag.eq(&other_name_trag));
        assert!(!other_name_trag.eq(&name_tag));
        assert!(name_tag.ne(&other_name_trag));
        assert!(other_name_trag.ne(&name_tag));

        // If id present, only consider id
        let full_tag = Tag::full("WGS 84", "EPSG", 7030);
        assert!(full_tag.eq(&full_tag));
        assert!(!full_tag.ne(&full_tag));
        assert!(!name_tag.eq(&full_tag));
        assert!(name_tag.ne(&full_tag));
        assert!(!full_tag.eq(&name_tag));
        assert!(full_tag.ne(&name_tag));

        let same_id_tag = Tag::full("WGS 1984", "EPSG", 7030);
        assert!(full_tag.eq(&same_id_tag));
        assert!(!full_tag.ne(&same_id_tag));
        assert!(same_id_tag.eq(&full_tag));
        assert!(!same_id_tag.ne(&full_tag));

        let other_id_tag = Tag::full("WGS 84", "geokit", 7030);
        assert!(!full_tag.eq(&other_id_tag));
        assert!(full_tag.ne(&other_id_tag));
        assert!(!other_id_tag.eq(&full_tag));
        assert!(other_id_tag.ne(&full_tag));
    }

    #[test]
    fn display() {
        let name_tag = Tag::name("WGS 84");
        assert_eq!(name_tag.to_string(), String::from("WGS 84"));

        let full_tag = Tag::full("WGS 84", "epsg", 7030);
        assert_eq!(full_tag.to_string(), "WGS 84 (epsg:7030)");
    }
}
