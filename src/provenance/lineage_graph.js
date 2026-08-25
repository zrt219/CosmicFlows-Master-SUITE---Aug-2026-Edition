/**
 * @file lineage_graph.js
 * @description W3C PROV-O JSON-LD Lineage Graph Engine.
 * Represents entities, activities, agents, and cryptographic data derivations.
 * 
 * @module provenance/lineage_graph
 */

export class ProvEntity {
  constructor(id, attributes = {}) {
    this.id = id;
    this.type = 'prov:Entity';
    this.attributes = attributes;
    this.wasGeneratedBy = null;
    this.wasDerivedFrom = [];
    this.sha256 = attributes.sha256 || null;
  }
}

export class ProvActivity {
  constructor(id, algorithm, parameters = {}) {
    this.id = id;
    this.type = 'prov:Activity';
    this.algorithm = algorithm;
    this.parameters = parameters;
    this.used = [];
    this.startTime = new Date().toISOString();
    this.endTime = null;
    this.wasAssociatedWith = null;
  }

  complete() {
    this.endTime = new Date().toISOString();
  }
}

export class LineageGraph {
  constructor(context = {}) {
    this.context = context;
    this.entities = new Map();
    this.activities = new Map();
    this.agents = new Map();
  }

  addEntity(id, attributes = {}) {
    const entity = new ProvEntity(id, attributes);
    this.entities.set(id, entity);
    return entity;
  }

  addActivity(id, algorithm, parameters = {}) {
    const activity = new ProvActivity(id, algorithm, parameters);
    this.activities.set(id, activity);
    return activity;
  }

  recordDerivation(targetEntityId, sourceEntityId, activityId) {
    const target = this.entities.get(targetEntityId);
    const source = this.entities.get(sourceEntityId);
    const activity = this.activities.get(activityId);

    if (target && source) {
      target.wasDerivedFrom.push(source.id);
    }
    if (activity && source) {
      activity.used.push(source.id);
    }
    if (target && activity) {
      target.wasGeneratedBy = activity.id;
    }
  }

  toJSONLD() {
    return {
      '@context': {
        'prov': 'http://www.w3.org/ns/prov#',
        'zrt': 'https://zrt.science/ontology/cosmicflows#'
      },
      'entities': Array.from(this.entities.values()),
      'activities': Array.from(this.activities.values()),
      'agents': Array.from(this.agents.values()),
      'generatedAt': new Date().toISOString()
    };
  }
}
