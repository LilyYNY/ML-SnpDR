"""Canonical identifiers shared by the Python stages of ML-SnpDR."""

from __future__ import annotations

from dataclasses import dataclass
import re


NETWORKS = {
    "string": ("String", "string"),
    "physicalppin": ("physicalPPIN", "physicalppin"),
    "chengf": ("chengF", "chengf"),
}
METHODS = {
    "louvain": ("Louvain", "louvain"),
    "wf": ("WF", "wf"),
}


def _token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(value).strip()).lower()


def _controlled(value: str, mapping: dict[str, tuple[str, str]], field: str) -> tuple[str, str]:
    key = _token(value)
    if key not in mapping:
        raise ValueError(f"Unknown {field}: {value}")
    return mapping[key]


@dataclass(frozen=True)
class ModuleIdentity:
    network: str
    method: str
    subtype: str
    module: str

    @property
    def uid(self) -> str:
        return make_module_uid(self.network, self.method, self.subtype, self.module)

    @property
    def legacy_module_id(self) -> str:
        network = _controlled(self.network, NETWORKS, "network")[0]
        method = _controlled(self.method, METHODS, "method")[0]
        return "|".join((network, method, self.subtype.upper(), self.module.upper()))


def make_module_uid(network: str, method: str, subtype: str, module: str) -> str:
    network_slug = _controlled(network, NETWORKS, "network")[1]
    method_slug = _controlled(method, METHODS, "method")[1]
    subtype_value = str(subtype).strip().upper()
    module_value = str(module).strip().upper()
    if not re.fullmatch(r"C[1-4]", subtype_value):
        raise ValueError("Subtype must be one of C1, C2, C3 or C4")
    if not re.fullmatch(r"M[0-9]+", module_value):
        raise ValueError("Module must use the form M<number>")
    return "__".join((network_slug, method_slug, subtype_value, module_value))


def parse_module_uid(module_uid: str) -> ModuleIdentity:
    fields = str(module_uid).split("__")
    if len(fields) != 4:
        raise ValueError(f"Invalid module_uid: {module_uid}")
    network = _controlled(fields[0], NETWORKS, "network")[0]
    method = _controlled(fields[1], METHODS, "method")[0]
    identity = ModuleIdentity(network, method, fields[2].upper(), fields[3].upper())
    if identity.uid != module_uid:
        raise ValueError(f"module_uid is not canonical: {module_uid}")
    return identity

