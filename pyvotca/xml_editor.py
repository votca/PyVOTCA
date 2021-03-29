"""Module to edit the options in the XML file."""

import xml.etree.ElementTree as ET
from typing import Any, Dict


def edit_xml(root: ET.Element, sections: Dict[str, Any], path: str = ".dftgwbse"):
    """Edit a XML object using sections."""
    sec = root.find(path)
    if sec is None:
        raise RuntimeError(f"Unkown Section: {path}")
    else:
        for key, val in sections.items():
            update_node(sec, key, val)


def update_node(root: ET.Element, key: str, value: Any):
    """Update nodes recursively."""
    sec = root.find(key)

    # insert new empty node
    if sec is None:
        elem = ET.Element(key)
        root.insert(0, elem)
        update_node(root, key, value)

    else:
        for node in root.findall(key):
            if not isinstance(value, dict):
                node.text = str(value)
            else:
                for k, v in value.items():
                    update_node(node, k, v)
