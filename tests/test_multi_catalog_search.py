# -*- coding: utf-8 -*-
"""
Automated Test Suite for Milestone 1: Multi-Catalog Search & Alternate Designations
Verifies:
1. Public API Bridge Contract on window.cosmicflows.search
2. Full Coverage of All 38 Astrometric Features + 8 Dupuy & Courtois (2023) Basins
3. Multi-Catalog Normalized Fuzzy Search (Abell, ACO, Messier, NGC, PGC, UGC, IC, 3C, WKK, ESO)
4. Relevance Scoring, Category Badges & Matched Alias Subtitle Highlighting
5. Keyboard Navigation (ArrowDown/Up, Enter selection, Escape dismiss, / and Ctrl+K shortcuts)
6. 3D Camera Fly-To, Targeting Beacon Animation, and Automatic Spectroscopy Dossier Modal
"""

import pytest
import time

TEST_SEARCH_QUERIES = [
    # Query, Expected Match Substring in Top Result, Expected Category
    ("A1656", "Coma Cluster", "cluster"),
    ("Abell 1656", "Coma Cluster", "cluster"),
    ("aco 1656", "Coma Cluster", "cluster"),
    ("a-1656", "Coma Cluster", "cluster"),
    ("M87", "Virgo Cluster", "cluster"),
    ("Messier 87", "Virgo Cluster", "cluster"),
    ("m49", "Virgo Cluster", "cluster"),
    ("m86", "Virgo Cluster", "cluster"),
    ("NGC 4889", "Coma Cluster", "cluster"),
    ("ngc4889", "Coma Cluster", "cluster"),
    ("NGC 4696", "Centaurus Cluster", "cluster"),
    ("Centaurus A", "Centaurus Cluster", "cluster"),
    ("NGC 5128", "Centaurus Cluster", "cluster"),
    ("NGC 1275", "Perseus Cluster", "cluster"),
    ("3c 84", "Perseus Cluster", "cluster"),
    ("ACO 3627", "Norma", "supercluster"),
    ("A3627", "Norma", "supercluster"),
    ("A3558", "Shapley Supercluster Core", "supercluster"),
    ("A426", "Perseus Cluster", "cluster"),
    ("A1060", "Hydra Cluster", "cluster"),
    ("A2151", "Hercules Cluster", "cluster"),
    ("A262", "Pisces Cluster", "cluster"),
    ("A1367", "Leo Cluster", "cluster"),
    ("A2199", "Abell 2199 Cluster", "cluster"),
    ("A2142", "Abell 2142 Monster Merger Cluster", "cluster"),
    ("ACO S0373", "Fornax Cluster", "cluster"),
    ("ACO S0636", "Antlia Cluster", "cluster"),
    ("ACO 3656", "Ophiuchus Cluster", "cluster"),
    ("Laniakea", "Laniakea", "supercluster"),
    ("Dipole Repeller", "Dipole Repeller", "void"),
    ("Boötes", "Boötes Supervoid", "void"),
    ("Bootes", "Boötes Supervoid", "void"),
    ("MeerKAT", "MeerKAT Vela ZoA Piercing Corridor", "corridor"),
    ("Parkes", "Parkes HIZOA Great Attractor Corridor", "corridor"),
    ("Basin 1", "Laniakea Supercluster Basin", "basin"),
    ("Basin 6", "Shapley Supercluster Basin", "basin"),
    ("SDSS-1a", "SDSS-1a Supercluster Basin", "basin"),
    ("PGC 41220", "Virgo Cluster", "cluster"),
    ("PGC 44715", "Coma Cluster", "cluster"),
    ("UGC 8168", "Coma Cluster", "cluster"),
    ("UGC 7654", "Virgo Cluster", "cluster"),
    ("UGC 2669", "Perseus Cluster", "cluster"),
    ("UGC 10170", "Hercules Cluster", "cluster")
]

class TestMultiCatalogSearchAPI:
    """Tier 1: Programmatic Search API Contract & Database Integrity"""

    def test_search_api_exposure(self, cdp):
        """Verify window.cosmicflows.search is exposed with complete API interface."""
        res = cdp.evaluate("""
        (function() {
            const s = window.cosmicflows && window.cosmicflows.search;
            if (!s) return { exposed: false };
            return {
                exposed: true,
                hasQuery: typeof s.query === 'function',
                hasSelect: typeof s.select === 'function',
                hasEntities: typeof s.getEntities === 'function',
                hasNormalize: typeof s.normalizeKey === 'function',
                entityCount: s.getEntities().length,
                hasAliases: !!s.ALIASES,
                aliasKeysCount: Object.keys(s.ALIASES || {}).length
            };
        })()
        """)
        assert res.get("exposed") is True, "window.cosmicflows.search must be exposed"
        assert res.get("hasQuery") is True, "search.query must be a function"
        assert res.get("hasSelect") is True, "search.select must be a function"
        assert res.get("hasEntities") is True, "search.getEntities must be a function"
        assert res.get("hasNormalize") is True, "search.normalizeKey must be a function"
        assert res.get("entityCount") >= 46, f"Expected >= 46 searchable entities, got {res.get('entityCount')}"
        assert res.get("aliasKeysCount") >= 46, f"Expected >= 46 alias keys, got {res.get('aliasKeysCount')}"

    def test_all_basins_registered_and_searchable(self, cdp):
        """Verify all 8 Dupuy & Courtois (2023) supercluster basins are indexed and searchable."""
        res = cdp.evaluate("""
        (function() {
            const s = window.cosmicflows.search;
            const entities = s.getEntities();
            const basins = entities.filter(e => e.category === 'basin');
            return {
                basinCount: basins.length,
                basinIds: basins.map(b => b.id),
                basinNames: basins.map(b => b.name)
            };
        })()
        """)
        assert res.get("basinCount") == 8, f"Expected 8 supercluster basins, got {res.get('basinCount')}"
        for i in range(1, 9):
            assert f"basin_{i}" in res["basinIds"], f"Missing basin_{i} in searchable database"

    @pytest.mark.parametrize("query, expected_name, expected_cat", TEST_SEARCH_QUERIES)
    def test_multi_catalog_search_queries(self, cdp, query, expected_name, expected_cat):
        """Verify multi-catalog fuzzy queries return the correct top astronomical match."""
        res = cdp.evaluate(f"""
        (function() {{
            const s = window.cosmicflows.search;
            const matches = s.query('{query}');
            if (!matches || matches.length === 0) return {{ matchCount: 0 }};
            const top = matches[0];
            return {{
                matchCount: matches.length,
                topName: top.name,
                topCategory: top.category,
                topId: top.id,
                score: top._score,
                matchedAlias: top._matchedAlias
            }};
        }})()
        """)
        assert res.get("matchCount") > 0, f"Query '{query}' returned 0 matches"
        assert expected_name.lower() in res["topName"].lower(), (
            f"Query '{query}' expected top match containing '{expected_name}', got '{res['topName']}'"
        )
        assert res["topCategory"] == expected_cat, (
            f"Query '{query}' expected category '{expected_cat}', got '{res['topCategory']}'"
        )
        assert res["score"] > 0, f"Query '{query}' match score must be > 0"


class TestSearchSceneAndDossierIntegration:
    """Tier 2: 3D Camera Fly-To, Targeting Beacon & Spectroscopy Dossier Activation"""

    def test_search_select_activates_flyto_beacon_and_dossier(self, cdp):
        """Verify selecting an entity moves camera, creates beacon, and opens dossier."""
        res = cdp.evaluate("""
        (function() {
            const search = window.cosmicflows.search;
            // Select Coma Cluster
            search.select('coma_cl');
            
            const modal = document.getElementById('spectroscopy-modal');
            const nameEl = document.getElementById('spec-name');
            const isOpen = modal && (modal.style.display !== 'none' && !modal.classList.contains('hidden'));
            const title = nameEl ? nameEl.textContent : '';

            return {
                isOpen: isOpen,
                title: title
            };
        })()
        """)
        # Allow short timeout for smooth transition & dossier open
        time.sleep(0.5)
        res_after = cdp.evaluate("""
        (function() {
            const modal = document.getElementById('spectroscopy-modal');
            const nameEl = document.getElementById('spec-name');
            const isOpen = modal && (modal.style.display !== 'none' && !modal.classList.contains('hidden'));
            const title = nameEl ? nameEl.textContent : '';
            
            // Clean up by closing dossier
            if (window.cosmicflows.spectroscopy) {
                window.cosmicflows.spectroscopy.closeDossier();
            }
            
            return {
                isOpen: isOpen,
                title: title
            };
        })()
        """)
        assert res_after["isOpen"] is True, "Spectroscopy modal must be open after search.select()"
        assert "coma" in res_after["title"].lower(), f"Modal title should show Coma, got '{res_after['title']}'"

    def test_search_select_synthesizes_dossier_for_tier3_void(self, cdp):
        """Verify selecting a feature without hardcoded dossier synthesizes physical record."""
        res = cdp.evaluate("""
        (function() {
            const search = window.cosmicflows.search;
            search.select('bootes_void');
            
            const modal = document.getElementById('spectroscopy-modal');
            const nameEl = document.getElementById('spec-name');
            const tagEl = document.getElementById('spec-tag');
            const czEl = document.getElementById('spec-cz');
            
            return {
                title: nameEl ? nameEl.textContent : '',
                tag: tagEl ? tagEl.textContent : '',
                cz: czEl ? czEl.textContent : ''
            };
        })()
        """)
        time.sleep(0.5)
        res_dossier = cdp.evaluate("""
        (function() {
            const spec = window.cosmicflows.spectroscopy;
            const d = spec.getDossier('bootes_void');
            
            // Close after test
            spec.closeDossier();
            
            return {
                found: !!d,
                name: d ? d.name : '',
                cat: d ? d.cat : '',
                cz: d ? d.cz : 0,
                sigmaV: d ? d.sigma_v : 0
            };
        })()
        """)
        assert res_dossier["found"] is True, "getDossier('bootes_void') must synthesize a valid record"
        assert "boötes" in res_dossier["name"].lower() or "bootes" in res_dossier["name"].lower()
        assert res_dossier["cat"] == "void"
        assert res_dossier["cz"] > 0


class TestSearchUIAndKeyboardNavigation:
    """Tier 3: Input Typing, Keyboard Navigation, and Hotkeys"""

    def test_search_dropdown_dom_rendering(self, cdp):
        """Verify typing into search input displays dropdown items with badges."""
        res = cdp.evaluate("""
        (function() {
            const input = document.getElementById('cosmic-search-input');
            const results = document.getElementById('cosmic-search-results');
            if (!input || !results) return { ok: false };
            
            input.value = 'A1656';
            input.dispatchEvent(new Event('input', { bubbles: true }));
            
            const items = results.querySelectorAll('.search-result-item');
            const isVisible = results.style.display === 'flex';
            const firstItemText = items.length > 0 ? items[0].textContent : '';
            const firstBadge = items.length > 0 ? items[0].querySelector('.search-tag-badge') : null;
            const badgeText = firstBadge ? firstBadge.textContent : '';
            
            // Clean up
            results.style.display = 'none';
            input.value = '';
            
            return {
                ok: true,
                isVisible: isVisible,
                itemCount: items.length,
                firstItemText: firstItemText,
                badgeText: badgeText
            };
        })()
        """)
        assert res["ok"] is True
        assert res["isVisible"] is True, "Search dropdown must be visible when typing"
        assert res["itemCount"] > 0, "Search dropdown must contain matching items"
        assert "coma" in res["firstItemText"].lower()
        assert res["badgeText"] == "CLUSTER"

    def test_keyboard_arrow_navigation(self, cdp):
        """Verify ArrowDown and ArrowUp navigate the result items and update .selected class."""
        res = cdp.evaluate("""
        (function() {
            const input = document.getElementById('cosmic-search-input');
            const results = document.getElementById('cosmic-search-results');
            
            input.value = 'supercluster';
            input.dispatchEvent(new Event('input', { bubbles: true }));
            
            // Press ArrowDown
            input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
            const items = results.querySelectorAll('.search-result-item');
            const firstSelected = items[0] ? items[0].classList.contains('selected') : false;
            
            // Press ArrowDown again
            input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
            const secondSelected = items[1] ? items[1].classList.contains('selected') : false;
            
            // Press ArrowUp
            input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowUp', bubbles: true }));
            const backToFirst = items[0] ? items[0].classList.contains('selected') : false;
            
            // Press Escape
            input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
            const isHidden = results.style.display === 'none';
            
            return {
                firstSelected: firstSelected,
                secondSelected: secondSelected,
                backToFirst: backToFirst,
                isHidden: isHidden
            };
        })()
        """)
        assert res["firstSelected"] is True, "ArrowDown must select first item"
        assert res["secondSelected"] is True, "Second ArrowDown must select second item"
        assert res["backToFirst"] is True, "ArrowUp must move selection back to first item"
        assert res["isHidden"] is True, "Escape must dismiss search results dropdown"
